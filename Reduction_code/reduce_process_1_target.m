function [finaldata]=reduce_process_1_target(objectname,options,raw_data_location,reduced_data_location,teff,logg,vsini)

%Initialise reduction state
finaldata=[];
reduction_performed=false;
skip_to_post_processing=false;

%REDUCTION MASTER FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Track starting directory
original_directory=pwd;
cleanup_object=onCleanup(@() cd(original_directory));

%Validate reduced-data location
if ~isfolder(reduced_data_location)
    error(['Reduced data location does not exist: ' reduced_data_location])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use exact target folder only.
%
%Do not search sibling folders such as:
%   objectname_old
%   objectname_2025
%   objectname_test
%
%Only files in reduced_data_location/objectname count as already reduced/rejected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reduced_target_directory=fullfile(reduced_data_location,objectname);
reduced_frames_directory=fullfile(reduced_target_directory,'reduced frames');
processing_directory=fullfile(reduced_frames_directory,'processing_files');
supplementary_directory=fullfile(reduced_target_directory,'supplementary_data');

%Create exact target folder if needed
if ~isfolder(reduced_target_directory)
    mkdir(reduced_target_directory);
end

%Create processing folder if needed
if ~isfolder(reduced_frames_directory), mkdir(reduced_frames_directory); end
if ~isfolder(processing_directory), mkdir(processing_directory); end
if ~isfolder(supplementary_directory), mkdir(supplementary_directory); end

%Migrate legacy target contents into the MEGARA 2.1 folder layout.
organise_existing_target_files(reduced_target_directory, ...
    reduced_frames_directory,processing_directory,supplementary_directory);

%Find existing reduced files, including rejected files, from the exact target folder only
past_reduction=list_accounted_reduction_files(reduced_frames_directory);

%Load star database
star_database=load_star_database_for_object(objectname,reduced_target_directory, ...
    supplementary_directory,reduced_data_location);
if isempty(star_database)
    disp(['No observations found for ' objectname '.'])
    finaldata=[];
    return
end

%Count observations
if ~isfield(star_database,'exptime')
    error(['star_database for ' objectname ' does not contain field exptime.'])
end

number_observations=length(star_database.exptime);
past_reduction_list=get_observation_numbers_from_reduction_files(past_reduction);
number_observations_reduced=length(past_reduction_list);

%Skip calibration objects
if strcmpi(objectname,'dark')
    disp(['Skipping reduction for ' objectname '.'])
    finaldata=[];
    return
end

if strcmpi(objectname,'thorium')
    disp(['Skipping reduction for ' objectname '.'])
    finaldata=[];
    return
end

%Report target status
fprintf('A total of %.0f observations have been found\n',number_observations)
fprintf('A total of %.0f observations have been previously reduced or rejected in exact target folder\n',number_observations_reduced)

%Skip reduction if requested
if options.skip_reduction
    disp(['Skipping reduction for ' objectname '.'])
else

    %Stop if there are no observations
    if number_observations==0
        disp(['Nothing to reduce for ' objectname ', please check objectname'])
        finaldata=[];
        return
    end

    %Check observation numbers in database
    database_observation_numbers=get_observation_numbers_from_database(star_database);

    %Decide which month/run folders need reduction
    if isempty(options.reduction_folders) && number_observations>length(past_reduction_list)

        if ~isfield(star_database,'months') || isempty(star_database.months)
            error(['star_database for ' objectname ' does not contain usable months information.'])
        end

        months=normalise_month_list(star_database.months);

    elseif ~isempty(options.reduction_folders)

        months=normalise_month_list(options.reduction_folders);

    else

        disp('No files to reduce, moving to post-reduction')
        skip_to_post_processing=true;

    end

    %Reduce required months
    if ~skip_to_post_processing

        for month_index=1:length(months)

            raw_month_directory=fullfile(raw_data_location,months{month_index});

            if ~isfolder(raw_month_directory)
                warning(['Raw data month/run folder does not exist: ' raw_month_directory])
                continue
            end

            cd(raw_month_directory)
            fprintf('Reducing Month/Run %.0f of %.0f, This is %s\n',month_index,length(months),months{month_index})

            %Load or create file_info
            if exist('./file_info.mat','file')==2
                loaded_file_information=load('./file_info.mat');

                if isfield(loaded_file_information,'info')
                    info=loaded_file_information.info;
                else
                    error(['file_info.mat in ' raw_month_directory ' does not contain variable info.'])
                end
            else
                info=display_directory_file_information_nooutput();
                save('file_info','info')
            end

            if numel(info.list)==0
                disp('No files to reduce, moving to next Month/Run')
                continue
            end

            %Find observation-number range in this raw folder
            first_file_name=info.list{1};
            last_file_name=info.list{end};

            start_observation_number=str2double(first_file_name(2:end));
            stop_observation_number=str2double(last_file_name(2:end));

            observations_in_month=database_observation_numbers>=start_observation_number & database_observation_numbers<=stop_observation_number;
            number_should_reduce=sum(observations_in_month);

            if isempty(past_reduction_list)
                number_already_reduced=0;
            else
                reductions_in_month=past_reduction_list>=start_observation_number & past_reduction_list<=stop_observation_number;
                number_already_reduced=sum(reductions_in_month);
            end

            if number_already_reduced~=number_should_reduce

                reduce_hercules_frames(objectname,past_reduction, ...
                    options.blue_data_chop_value,options.skip_manual_file_check)
                reduction_performed=true;

                %Move newly reduced files from raw month folder into exact reduced target folder
                move_files_matching_pattern(raw_month_directory,'J*.mat',reduced_frames_directory);

                %Move processing files into exact processing folder
                move_files_matching_pattern(reduced_frames_directory,'J*_prc.mat',processing_directory);

                %Refresh accounted reductions from exact target folder only
                past_reduction=list_accounted_reduction_files(reduced_frames_directory);
                past_reduction_list=get_observation_numbers_from_reduction_files(past_reduction);

            else

                fprintf('All %.0f observations in %s are already reduced or rejected in exact target folder\n',number_should_reduce,months{month_index})

            end
        end
    end

    %Move any processing files after reduction
    cd(reduced_target_directory)
    move_files_matching_pattern(reduced_frames_directory,'J*_prc.mat',processing_directory);

end

%Move to exact reduced target folder
cd(reduced_target_directory)

%Move processing files before post-reduction
move_files_matching_pattern(reduced_frames_directory,'J*_prc.mat',processing_directory);

%Find final active reduced files only. Rejected files must not be re-entered into post-reduction.
final_files=list_active_reduction_files(reduced_frames_directory);
final_data_file_name=get_final_data_file_name(reduced_target_directory,objectname,options);
process_data_flag=false;

%Run processing and post-reduction
if ~options.skip_post_reduction

    if isempty(final_files)

        finaldata=[];
        process_data_flag=false;

    else

        finaldata=load_final_data_file(final_data_file_name);

        if ~options.overwrite_full_data && ~reduction_performed && ~options.calculate_lsd_profile && ~options.radial_velocity_measurement && ~options.FAMIAS_output && ~isempty(finaldata)

            disp('Star previously processed')

            if ~options.suppress_figures
                plot_final_spectrum(finaldata,objectname)
            end

        else

            process_data_flag=should_process_final_data(finaldata,options,reduced_target_directory,objectname);

            if process_data_flag
                finaldata=process_1_target(objectname,options,reduced_data_location,teff,logg,vsini);
            end

        end
    end

    if ~options.overwrite_post_reduction || ~process_data_flag
        disp('Skipping post-reduction')
    end

else

    finaldata=load_final_data_file(final_data_file_name);
    disp('Skipping processing and post-reduction')

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function move_files_matching_pattern(source_directory,file_pattern,destination_directory)

%Move matching files
if ~isfolder(destination_directory)
    mkdir(destination_directory)
end

matching_files=dir(fullfile(source_directory,file_pattern));

for matching_file_index=1:numel(matching_files)

    if ~matching_files(matching_file_index).isdir

        source_file_name=fullfile(source_directory,matching_files(matching_file_index).name);
        [move_status,move_message]=movefile(source_file_name,destination_directory,'f');

        if ~move_status
            warning(['Could not move file ' source_file_name ': ' move_message])
        end

    end
end

end

function organise_existing_target_files(target_directory,reduced_frames_directory, ...
    processing_directory,supplementary_directory)

legacy_processing_directory=fullfile(target_directory,'processing_files');
move_files_matching_pattern(target_directory,'J*_prc.mat',processing_directory);
move_files_matching_pattern(reduced_frames_directory,'J*_prc.mat',processing_directory);
if isfolder(legacy_processing_directory)
    move_directory_contents(legacy_processing_directory,processing_directory);
end
move_files_matching_pattern(target_directory,'J*.mat',reduced_frames_directory);
move_directory_contents(fullfile(target_directory,'removed_reduced_files'), ...
    fullfile(reduced_frames_directory,'removed_reduced_files'));
move_directory_contents(fullfile(target_directory,'removed_files'), ...
    fullfile(reduced_frames_directory,'removed_files'));

information_patterns={'file_info.mat','star_database_*.mat','datacf_*.mat', ...
    'post_reduction_*.mat'};
for pattern_index=1:numel(information_patterns)
    move_files_matching_pattern(target_directory,information_patterns{pattern_index}, ...
        supplementary_directory);
end

move_directory_contents(fullfile(target_directory,'automatic_continuum_results'), ...
    fullfile(supplementary_directory,'automatic_continuum_results'));
move_directory_contents(fullfile(target_directory,'chunked_post_reduction_files'), ...
    fullfile(supplementary_directory,'chunked_post_reduction_files'));

end

function move_directory_contents(source_directory,destination_directory)

if ~isfolder(source_directory), return; end
if ~isfolder(destination_directory), mkdir(destination_directory); end
directory_contents=dir(source_directory);
for item_index=1:numel(directory_contents)
    item_name=directory_contents(item_index).name;
    if any(strcmp(item_name,{'.','..'})), continue; end
    movefile(fullfile(source_directory,item_name),destination_directory,'f');
end

end

function star_database=load_star_database_for_object(objectname,reduced_target_directory,supplementary_directory,reduced_data_location)

%Initialise database
star_database=[];

%Build possible database names.
%This searches the exact target folder, the observation-summary folder, and the reduced-data root.
%It does not search sibling reduced folders such as objectname_old or objectname_2025.
observation_summary_directory=fullfile(reduced_data_location,'Observation Summaries');

candidate_database_files={ ...
    fullfile(supplementary_directory,['star_database_' objectname '.mat']), ...
    fullfile(reduced_target_directory,['star_database_' objectname '.mat']), ...
    fullfile(observation_summary_directory,['star_database_' objectname '.mat']), ...
    fullfile(reduced_data_location,['star_database_' objectname '.mat']) ...
};

if numel(objectname)>2
    candidate_database_files{end+1}=fullfile(supplementary_directory,['star_database_HD' objectname(3:end) '.mat']);
    candidate_database_files{end+1}=fullfile(reduced_target_directory,['star_database_HD' objectname(3:end) '.mat']);
    candidate_database_files{end+1}=fullfile(observation_summary_directory,['star_database_HD' objectname(3:end) '.mat']);
    candidate_database_files{end+1}=fullfile(reduced_data_location,['star_database_HD' objectname(3:end) '.mat']);
end

%Load first matching database
for candidate_file_index=1:numel(candidate_database_files)

    if isfile(candidate_database_files{candidate_file_index})
        loaded_database=load(candidate_database_files{candidate_file_index});
        star_database=get_database_variable_from_loaded_file(loaded_database,objectname,candidate_database_files{candidate_file_index});
        return
    end

end

end

function star_database=get_database_variable_from_loaded_file(loaded_database,objectname,database_file_name)

%Extract database variable
expected_variable_name=['star_database_' objectname];

if isfield(loaded_database,expected_variable_name)

    star_database=loaded_database.(expected_variable_name);

elseif numel(objectname)>2 && isfield(loaded_database,['star_database_HD' objectname(3:end)])

    star_database=loaded_database.(['star_database_HD' objectname(3:end)]);

elseif isfield(loaded_database,'star_database')

    star_database=loaded_database.star_database;

else

    loaded_variable_names=fieldnames(loaded_database);
    star_database=loaded_database.(loaded_variable_names{1});
    warning(['Could not find expected star database variable in ' database_file_name '. Using first variable in file instead.'])

end

end

function observation_numbers=get_observation_numbers_from_reduction_files(reduction_files)

%Extract reduced observation numbers
observation_numbers=[];

for reduction_file_index=1:numel(reduction_files)

    reduction_file_name=reduction_files(reduction_file_index).name;
    observation_match=regexp(reduction_file_name,'^J(\d+)\.mat$','tokens','once');

    if isempty(observation_match)
        continue
    end

    observation_numbers(end+1)=str2double(observation_match{1}); %#ok<AGROW>

end

observation_numbers=unique(observation_numbers);

end

function database_observation_numbers=get_observation_numbers_from_database(star_database)

%Extract database observation numbers
if ~isfield(star_database,'filename')
    error('star_database does not contain field filename.')
end

filename_list=star_database.filename;
database_observation_numbers=nan(numel(filename_list),1);

for filename_index=1:numel(filename_list)

    single_filename=filename_list{filename_index};
    observation_match=regexp(single_filename,'^J(\d+)\.fit$','tokens','once');

    if isempty(observation_match)
        database_observation_numbers(filename_index)=NaN;
    else
        database_observation_numbers(filename_index)=str2double(observation_match{1});
    end

end

database_observation_numbers=database_observation_numbers(isfinite(database_observation_numbers));
database_observation_numbers=unique(database_observation_numbers);

end

function months=normalise_month_list(months)

%Normalise month list
if ischar(months)
    months={months};
end

if isstring(months)
    months=cellstr(months);
end

if ~iscellstr(months)
    error('months must be a cell array of character vectors, string array, character vector, or empty.')
end

months=unique(months,'stable');

end

function final_data_file_name=get_final_data_file_name(reduced_target_directory,objectname,options)

%Choose final-data filename
if options.order_merge
    final_data_file_name=fullfile(reduced_target_directory,['final_data_' objectname '.mat']);
else
    final_data_file_name=fullfile(reduced_target_directory,['final_data_unmerged_' objectname '.mat']);
end

end

function finaldata=load_final_data_file(final_data_file_name)

%Load final data
finaldata=[];

if isfile(final_data_file_name)

    loaded_final_data=load(final_data_file_name);

    if isfield(loaded_final_data,'finaldata')
        finaldata=loaded_final_data.finaldata;
    else
        loaded_variable_names=fieldnames(loaded_final_data);
        finaldata=loaded_final_data.(loaded_variable_names{1});
        warning(['Could not find variable finaldata in ' final_data_file_name '. Using first variable in file instead.'])
    end

end

end

function process_data_flag=should_process_final_data(finaldata,options,reduced_target_directory,objectname)

%Decide whether to process
regularised_lsd_profile_file=fullfile(reduced_target_directory,['regularised_lsd_profile_' objectname '.mat']);

if isempty(finaldata)

    process_data_flag=true;

elseif options.overwrite_post_reduction || options.overwrite_full_data

    process_data_flag=true;

elseif options.calculate_lsd_profile && ~isfile(regularised_lsd_profile_file)

    process_data_flag=true;

elseif options.radial_velocity_measurement && (~isfield(finaldata,'radial_velocity') || all(isnan(finaldata.radial_velocity(:))))

    process_data_flag=true;

elseif isfield(finaldata,'systemic_velocity_direct') && isnumeric(finaldata.systemic_velocity_direct) && isscalar(finaldata.systemic_velocity_direct) && isnan(finaldata.systemic_velocity_direct)

    process_data_flag=true;

else

    process_data_flag=false;

end

end

function plot_final_spectrum(finaldata,objectname)

%Plot final spectrum
if isfield(finaldata,'fullwave') && isfield(finaldata,'fullint')

    figure
    plot(finaldata.fullwave,finaldata.fullint)
    title(['Full Spectra ' get_object_display_name(objectname)])

else

    warning(['Cannot plot final spectrum for ' objectname ' because finaldata.fullwave or finaldata.fullint is missing.'])

end

end

function object_display_name=get_object_display_name(objectname)

%Format object label
if numel(objectname)>=4
    object_display_name=objectname(4:end);
else
    object_display_name=objectname;
end

end

function active_reduction_files=list_active_reduction_files(reduced_frames_directory)

%List active reduced files only from exact target folder
active_reduction_files=dir(fullfile(reduced_frames_directory,'J*.mat'));
keep_file=false(numel(active_reduction_files),1);

for file_index=1:numel(active_reduction_files)

    keep_file(file_index)=~active_reduction_files(file_index).isdir && ~isempty(regexp(active_reduction_files(file_index).name,'^J\d+\.mat$','once'));

end

active_reduction_files=active_reduction_files(keep_file);

end

function accounted_reduction_files=list_accounted_reduction_files(reduced_frames_directory)

%List reduced files that should count as already handled.
%This includes active reduced files and reduced files deliberately rejected
%during post-reduction, so rejected observations are not regenerated every run.
%
%Important:
%Only the exact target folder and its own rejected-file folders are searched.
%Sibling folders such as objectname_old or objectname_2025 are deliberately ignored.
search_directories={ ...
    reduced_frames_directory, ...
    fullfile(reduced_frames_directory,'removed_files'), ...
    fullfile(reduced_frames_directory,'removed_reduced_files') ...
};

accounted_reduction_files=struct('name',{},'folder',{},'date',{},'bytes',{},'isdir',{},'datenum',{});

for directory_index=1:numel(search_directories)

    this_directory=search_directories{directory_index};

    if ~isfolder(this_directory)
        continue
    end

    this_files=dir(fullfile(this_directory,'J*.mat'));

    for file_index=1:numel(this_files)

        if this_files(file_index).isdir
            continue
        end

        %Keep only true reduced files, not processing files
        if isempty(regexp(this_files(file_index).name,'^J\d+\.mat$','once'))
            continue
        end

        accounted_reduction_files(end+1)=this_files(file_index); %#ok<AGROW>

    end
end

%Remove duplicate names, keeping first occurrence
if ~isempty(accounted_reduction_files)

    [~,unique_index]=unique({accounted_reduction_files.name},'stable');
    accounted_reduction_files=accounted_reduction_files(unique_index);

end

end
