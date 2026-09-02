function []=reduce_hercules_frames(objectname,past_red,blue_data_chop_value,skip_manual_file_check)

if nargin<4 || isempty(skip_manual_file_check)
    skip_manual_file_check=false;
end

%Track raw month/run directory
reduction_directory=pwd;
cleanup_object=onCleanup(@() cd(reduction_directory));

%Define cache files using absolute paths
flat_cache_file=fullfile(reduction_directory,['flat_info_blue_' num2str(blue_data_chop_value) '.mat']);
legacy_flat_cache_file=fullfile(reduction_directory, ...
    ['flat_info_blue_' num2str(blue_data_chop_value) '_Megara_1.6.mat']);
thar_cache_file=fullfile(reduction_directory,['ThAr_info_blue_' num2str(blue_data_chop_value) '.mat']);

%1-program to identify and remove with incomplete/mislabelled data. This runs automatically unless badfiles.mat exisits in the diectory
cd(reduction_directory)
[change_flag,counter_good]=file_check_routine( ...
    objectname,skip_manual_file_check);

%Return to raw month/run directory in case file checking changed folder
cd(reduction_directory)

%Load or rebuild file information
if change_flag==1
    fprintf(['Rejected files were moved; refreshing the FITS header ' ...
        'information for this Month/Run\n'])
    info=display_directory_file_information_nooutput();
    save(fullfile(reduction_directory,'file_info.mat'),'info')
else
    file_info_name=fullfile(reduction_directory,'file_info.mat');
    if isfile(file_info_name)
        loaded_file_info=load(file_info_name,'info');
        if isfield(loaded_file_info,'info')
            info=loaded_file_info.info;
        else
            error('file_info.mat exists but does not contain variable info.');
        end
    else
        info=display_directory_file_information_nooutput();
        save(file_info_name,'info')
    end
end

%Report matching frames without implying that absent frames failed quality checks
fprintf('A total of %.0f matching stellar frames are ready for reduction\n',counter_good)

%Stop if no matching stellar files need processing
if counter_good==0
    disp('No new matching stellar frames require reduction in this Month/Run')
    return
end

%A stellar extraction cannot be built without a compatible flat field.
number_white_lamps=0;
for file_index=1:numel(info.list)
    file_identifier=info.list{file_index};
    if isfield(info,file_identifier) && ...
            isfield(info.(file_identifier),'HERCEXPT') && ...
            strcmp(info.(file_identifier).HERCEXPT,'White L')
        number_white_lamps=number_white_lamps+1;
    end
end
if number_white_lamps==0 && ~isfile(flat_cache_file) && ...
        ~isfile(legacy_flat_cache_file)
    warning('MEGARA:MissingFlatCalibration', ...
        ['Skipping %.0f stellar frames in %s because this run has no ' ...
        'White L flat fields and no compatible current or legacy flat cache.'], ...
        counter_good,reduction_directory)
    return
end

%Set reduction parameters
bias_frame=[];
plotting=0; % 0 for no, 1 for yes flat field cosmic ray filtering plotting (default=30)
flatdata_medfilt=1; % amount of median smoothing on flat-field data (former default=30, now set to 1 to remove chip artefacts)
backdata_medfilt=50; % amount of mean smoothing on stellar and flat-field background data (default=50)
cosmic_plotting=0; % whether to plot the cosmic ray filtering process for stellar images (default=0)
number_bad_thorium=0; %number of thorium files with less than 900 lines found

%2 process the flat-fields by summing all flats, order-tracing and background fitting
cd(reduction_directory)
[allorders,flatdata,backdata,summed_flat,backfit,save_flag_flat_field]=flat_processing(info,bias_frame,blue_data_chop_value,backdata_medfilt,plotting);

%Save flat-field information if rebuilt
if save_flag_flat_field==1
    cd(reduction_directory)
    save(flat_cache_file,'allorders','flatdata','backdata','summed_flat','backfit')
end

%3- The following will process all the thar images taken that night
cd(reduction_directory)
[thorium_julian_dates,thorium_solution_metadata,thorium_fitted_wavelengths,number_bad_thorium,save_flag_thorium,sum_info_rebuild]=thar_processing(info,allorders,bias_frame,flatdata,blue_data_chop_value,number_bad_thorium);

%Save ThAr information if rebuilt
if save_flag_thorium==1
    cd(reduction_directory)
    save_thorium_information(thar_cache_file,thorium_julian_dates,thorium_solution_metadata,thorium_fitted_wavelengths)
end

%Rebuild file information if ThAr processing changed files
if sum_info_rebuild>=1
    cd(reduction_directory)
    delete_file_if_present(fullfile(reduction_directory,'file_information.txt'));
    delete_file_if_present(fullfile(reduction_directory,'file_info.mat'));
    info=display_directory_file_information_nooutput();
    save(fullfile(reduction_directory,'file_info.mat'),'info')
end

%4- Process the stellar images one by one using the two thoriums either side of the image (unless there aren't any!).
%include alias names

%Return to raw month/run directory before stellar processing
cd(reduction_directory)

%Build alternative star-name list
alt_starnames=get_alternative_star_names(objectname);

%Process each available star name
for alternative_name_index=1:numel(alt_starnames)

    starname=alt_starnames{alternative_name_index};

    %only process if have files of this name
    process_flag=false;

    for file_info_index=1:numel(info.list)
        [type,star]=identify_the_file(info,file_info_index); %#ok<ASGLU>
        if strcmp(star,starname)
            process_flag=true;
        end
    end

    if process_flag==1
        cd(reduction_directory)
        stellar_processing(info,allorders,flatdata,thorium_julian_dates,thorium_fitted_wavelengths,thorium_solution_metadata,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,starname,past_red);
    end

end

%notify the user that the run has finished
disp(' ')
disp('Run reduced')

if number_bad_thorium~=0
    warning(' %.0f thoriums may not have reduced correctly!\n',number_bad_thorium)
end

disp(' ')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_thorium_information(thar_cache_file,thorium_julian_dates,thorium_solution_metadata,thorium_fitted_wavelengths)

%Save legacy ThAr variable names
jd_ths=thorium_julian_dates;
m=thorium_solution_metadata;
fwave=thorium_fitted_wavelengths;

%Mark cache as MEGARA 2.0-compatible
megara_cache_version='MEGARA_2.0';

if ~isstruct(fwave)
    error('thorium_fitted_wavelengths must be a structure before saving ThAr cache.')
end

fwave.megara_cache_version=megara_cache_version;

%Save ThAr cache using explicit path
save(thar_cache_file,'jd_ths','m','fwave','megara_cache_version')

end

function delete_file_if_present(file_name)

%Delete file if present
if isfile(file_name)
    delete(file_name)
end

end

function alt_starnames=get_alternative_star_names(objectname)

%Initialise alternative names
alt_starnames={objectname};

%Load known aliases
if ~isfile('Mus_known_alias_names.mat')
    return
end

loaded_aliases=load('Mus_known_alias_names.mat');

if ~isfield(loaded_aliases,'wildc')
    return
end

wildc=loaded_aliases.wildc;

%Find aliases
alias_rows=find(any(strcmpi(wildc,objectname),2));

if ~isempty(alias_rows)
    for alias_row_index=1:numel(alias_rows)

        alias_row=alias_rows(alias_row_index);

        for alias_column_index=1:size(wildc,2)

            alias_name=wildc{alias_row,alias_column_index};

            if isempty(alias_name)
                continue
            end

            if ~any(strcmpi(alt_starnames,alias_name))
                alt_starnames{end+1}=alias_name; %#ok<AGROW>
            end

        end
    end
end

end
