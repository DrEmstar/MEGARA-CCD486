function write_ascii_output(objectname,final_data,barycentric_correction_flag,radial_velocity_measurement,calculate_lsd_profile)

%Store starting directory
starting_directory=pwd;
cleanup_directory=onCleanup(@() cd(starting_directory));

%Get target names
target_names=get_target_names_from_objectname(objectname);

%Write ASCII output for each target
for target_index=1:numel(target_names)
    target_name=target_names{target_index};
    final_data_field_name=matlab.lang.makeValidName(target_name);
    if ~isfield(final_data,final_data_field_name) || isempty(final_data.(final_data_field_name)), continue; end
    target_directory=fullfile(starting_directory,target_name);
    if ~isfolder(target_directory), warning(['Reduced target directory does not exist: ' target_directory]); continue; end
    cd(target_directory)
    target_data=final_data.(final_data_field_name);

    %Write regularised LSD profiles in dat format
    if calculate_lsd_profile==1
        write_regularised_lsd_ascii_file(target_name);
    end

    %Write order-by-order or merged spectra
    target_data_field_names=fieldnames(target_data);
    if any(startsWith(target_data_field_names,'int'))
        write_order_by_order_ascii_files(target_name,target_data,barycentric_correction_flag);
    elseif any(strcmp(target_data_field_names,'fullwave')) && any(strcmp(target_data_field_names,'fullint'))
        write_merged_ascii_files(target_name,target_data,barycentric_correction_flag,radial_velocity_measurement);
    else
        warning(['No recognised reduced spectra fields found for ' target_name '.'])
    end
    cd(starting_directory)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function target_names=get_target_names_from_objectname(objectname)

%Read target names
if isstruct(objectname)
    object_field_names=fieldnames(objectname);
    target_names=cell(numel(object_field_names),1);
    for object_index=1:numel(object_field_names)
        target_names{object_index}=objectname.(object_field_names{object_index});
    end
elseif ischar(objectname) || isstring(objectname)
    target_names=cellstr(string(objectname));
elseif isempty(objectname)
    target_names={};
else
    error('objectname must be a structure, character vector, string scalar, or empty.')
end

%Normalise target names
for target_index=1:numel(target_names)
    if isstring(target_names{target_index}), target_names{target_index}=char(target_names{target_index}); end
    if ~ischar(target_names{target_index}), error('Each object name must be a character vector or string scalar.'); end
end

end

function write_regularised_lsd_ascii_file(target_name)

%Write regularised LSD profiles in dat format
regularised_lsd_file_name=['regularised_lsd_profile_' target_name '.mat'];
if ~isfile(regularised_lsd_file_name), warning(['Regularised LSD profile file not found: ' regularised_lsd_file_name]); return; end

loaded_lsd_file=load(regularised_lsd_file_name);
if ~isfield(loaded_lsd_file,'regularised_lsd_profile'), warning(['Regularised LSD profile file does not contain regularised_lsd_profile: ' regularised_lsd_file_name]); return; end

regularised_lsd_profile=loaded_lsd_file.regularised_lsd_profile;
if ~isfield(regularised_lsd_profile,'velocity') || ~isfield(regularised_lsd_profile,'intensity') || ~isfield(regularised_lsd_profile,'jd'), warning(['regularised_lsd_profile is missing velocity, intensity, or jd in ' regularised_lsd_file_name]); return; end

[velocity,intensity]=standardise_lsd_profile_arrays(regularised_lsd_profile.velocity,regularised_lsd_profile.intensity,regularised_lsd_file_name);
jd=regularised_lsd_profile.jd(:);

if numel(jd)~=size(intensity,1)
    warning(['Number of JD values does not match number of regularised LSD profiles in ' regularised_lsd_file_name '. Writing available intensity profiles.'])
end

ascii_regularised_lsd_file=zeros(numel(velocity),size(intensity,1)+1);
ascii_regularised_lsd_file(:,1)=velocity(:);

for observation_index=1:size(intensity,1)
    ascii_regularised_lsd_file(:,observation_index+1)=intensity(observation_index,:)';
end

ascii_regularised_lsd_filename=['regularised_lsd_profile_' target_name '.dat'];
save(ascii_regularised_lsd_filename,'ascii_regularised_lsd_file','-ascii')

end

function [velocity,intensity]=standardise_lsd_profile_arrays(velocity,intensity,label_text)

%Force velocity vector into row orientation
velocity=velocity(:).';

%Check intensity matrix orientation
if size(intensity,2)==numel(velocity)
    %Already spectra x velocity
elseif size(intensity,1)==numel(velocity)
    intensity=intensity.';
else
    error([label_text ' has incompatible velocity/intensity sizes. velocity has ' num2str(numel(velocity)) ' points; intensity is ' num2str(size(intensity,1)) ' x ' num2str(size(intensity,2)) '.'])
end

%Sort by velocity if needed
[velocity,sort_index]=sort(velocity);
intensity=intensity(:,sort_index);

end

function write_order_by_order_ascii_files(target_name,target_data,barycentric_correction_flag)

%Find reduced files
file_list=list_reduced_spectrum_files();
if isempty(file_list), warning(['No reduced J*.mat files found for ' target_name]); return; end

%Find available orders
order_numbers=get_available_order_numbers(target_data);
if isempty(order_numbers), warning(['No order-by-order data found for ' target_name]); return; end
order_numbers=sort(order_numbers,'descend');

%Write one ASCII file per observation
summary_file=cell(numel(file_list),9);
for observation_index=1:numel(file_list)
    ascii_file=[];
    reduced_file_root=file_list(observation_index).name(1:end-4);
    reduced_frame=load_primary_variable_from_mat_file(file_list(observation_index).name);
    for order_index=1:numel(order_numbers)
        order_number=order_numbers(order_index);
        wave_field_name=['wave' num2str(order_number)];
        intensity_field_name=['int' num2str(order_number)];
        weights_field_name=['weights' num2str(order_number)];
        if ~isfield(target_data,wave_field_name) || ~isfield(target_data,intensity_field_name) || ~isfield(target_data,weights_field_name), continue; end
        wavelength_axis=target_data.(wave_field_name)(:);
        intensity_axis=target_data.(intensity_field_name)(observation_index,:)';
        weight_axis=target_data.(weights_field_name)(observation_index,:)';
        order_axis=order_number.*ones(numel(wavelength_axis),1);
        ascii_file_part=cat(2,wavelength_axis,intensity_axis,weight_axis,order_axis);
        ascii_file=cat(1,ascii_file,ascii_file_part);
    end
    summary_file(observation_index,:)=build_summary_row(reduced_file_root,reduced_frame,barycentric_correction_flag);
    ascii_filename=[reduced_file_root '_reduced.dat'];
    file_identifier=fopen(ascii_filename,'wt');
    if file_identifier==-1, error(['Could not open ASCII output file: ' ascii_filename]); end
    cleanup_file=onCleanup(@() fclose(file_identifier));
    fprintf(file_identifier,'Wavelength Intensity Weight Order \n');
    for row_index=1:size(ascii_file,1)
        fprintf(file_identifier,'%.2f %.12f %.12f %.0f\n',ascii_file(row_index,:));
    end
    clear cleanup_file
end

%Write summary file
summary_filename=['summary_file_' target_name '.txt'];
summary_table=cell2table(summary_file,'VariableNames',{'Filename','JD','Date','Exp_Start','Exp_Time','FWMT','BCorr','BCorr_App','Signal_to_Noise'});
writetable(summary_table,summary_filename)

end

function write_merged_ascii_files(target_name,target_data,barycentric_correction_flag,radial_velocity_measurement)

%Find reduced files
file_list=list_reduced_spectrum_files();
if isempty(file_list), warning(['No reduced J*.mat files found for ' target_name]); return; end

%Write one ASCII file per observation
if radial_velocity_measurement==0
    summary_file=cell(numel(file_list),9);
else
    summary_file=cell(numel(file_list),14);
end

for observation_index=1:numel(file_list)
    reduced_file_root=file_list(observation_index).name(1:end-4);
    reduced_frame=load_primary_variable_from_mat_file(file_list(observation_index).name);
    ascii_file=cat(2,target_data.fullwave(:),target_data.fullint(observation_index,:)');
    summary_file(observation_index,1:9)=build_summary_row(reduced_file_root,reduced_frame,barycentric_correction_flag);
    ascii_filename=[reduced_file_root '_reduced.dat'];
    file_identifier=fopen(ascii_filename,'wt');
    if file_identifier==-1, error(['Could not open ASCII output file: ' ascii_filename]); end
    cleanup_file=onCleanup(@() fclose(file_identifier));
    fprintf(file_identifier,'Wavelength Intensity \n');
    for row_index=1:size(ascii_file,1)
        fprintf(file_identifier,'%.2f %.12f\n',ascii_file(row_index,:));
    end
    clear cleanup_file
    if radial_velocity_measurement==1
        summary_file{observation_index,10}=get_scalar_or_nan(target_data,'systemic_velocity_gauss');
        summary_file{observation_index,11}=get_indexed_value_or_nan(target_data,'radial_velocity',observation_index);
        summary_file{observation_index,12}=get_indexed_or_scalar_value_or_nan(target_data,'radial_velocity_error',observation_index);
        summary_file{observation_index,13}=get_indexed_value_or_nan(target_data,'vsini',observation_index);
        summary_file{observation_index,14}=get_indexed_or_scalar_value_or_nan(target_data,'vsini_error',observation_index);
        summary_file{observation_index,15}=get_vbroad_value(target_data,observation_index);
        summary_file{observation_index,16}=get_vbroad_error(target_data,observation_index);
    end
end

%Write summary file
summary_filename=['summary_file_' target_name '.txt'];
if radial_velocity_measurement==0
    summary_table=cell2table(summary_file,'VariableNames',{'Filename','JD','Date','Exp_Start','Exp_Time','FWMT','BCorr','BCorr_App','Signal_to_Noise'});
else
    summary_table=cell2table(summary_file,'VariableNames',{'Filename','JD','Date','Exp_Start','Exp_Time','FWMT','BCorr','BCorr_App','Signal_to_Noise','Syst_Vel','RV','RV_err','Vsini','Vsini_err','Vbroad','Vbroad_err'});
end
writetable(summary_table,summary_filename)

end

function value=get_vbroad_value(target_data,observation_index)

%Use the descriptive alias, with compatibility for older final_data files.
if isfield(target_data,'vbroad')
    value=get_indexed_value_or_nan(target_data,'vbroad',observation_index);
else
    value=get_indexed_value_or_nan(target_data,'vsini',observation_index);
end

end

function value=get_vbroad_error(target_data,observation_index)

if isfield(target_data,'vbroad_error')
    value=get_indexed_or_scalar_value_or_nan(target_data,'vbroad_error',observation_index);
else
    value=get_indexed_or_scalar_value_or_nan(target_data,'vsini_error',observation_index);
end

end

function summary_row=build_summary_row(reduced_file_root,reduced_frame,barycentric_correction_flag)

%Initialise summary row
summary_row=cell(1,9);
summary_row{1}=reduced_file_root;
summary_row{2}=get_structure_field_or_default(reduced_frame,'jd',NaN);

%Read header metadata if FITS file is available
fits_file_name=[reduced_file_root '.fit'];
if isfile(fits_file_name)
    header=get_header(fits_file_name);
    headers=get_required_headers_from_header(header);
    summary_row{3}=get_structure_field_or_default(headers,'DATE_OBS',get_structure_field_or_default(reduced_frame,'expdate',''));
    summary_row{4}=get_structure_field_or_default(headers,'REC_STRT','');
    summary_row{5}=get_structure_field_or_default(headers,'EXPTIME',get_structure_field_or_default(reduced_frame,'exptime',NaN));
    summary_row{6}=get_structure_field_or_default(headers,'HERCFWMT',NaN);
else
    summary_row{3}=get_structure_field_or_default(reduced_frame,'expdate','');
    summary_row{4}=get_structure_field_or_default(reduced_frame,'midtime_utc','');
    summary_row{5}=get_structure_field_or_default(reduced_frame,'exptime',NaN);
    summary_row{6}=NaN;
end

%Read stored reduction metadata
summary_row{7}=get_structure_field_or_default(reduced_frame,'bcorr',NaN);
summary_row{8}=barycentric_correction_flag;
summary_row{9}=get_structure_field_or_default(reduced_frame,'signal_to_noise',NaN);

end

function order_numbers=get_available_order_numbers(target_data)

%Find order numbers from wave fields
target_data_field_names=fieldnames(target_data);
order_numbers=[];

for field_index=1:numel(target_data_field_names)
    token=regexp(target_data_field_names{field_index},'^wave(\d+)$','tokens','once');
    if isempty(token), continue; end
    order_number=str2double(token{1});
    intensity_field_name=['int' num2str(order_number)];
    weights_field_name=['weights' num2str(order_number)];
    if isfield(target_data,intensity_field_name) && isfield(target_data,weights_field_name)
        order_numbers(end+1)=order_number;
    end
end

order_numbers=unique(order_numbers);

end

function reduced_files=list_reduced_spectrum_files()

%List reduced stellar files only
all_files=dir('J*.mat');
keep_file=false(numel(all_files),1);

for file_index=1:numel(all_files)
    keep_file(file_index)=~isempty(regexp(all_files(file_index).name,'^J\d+\.mat$','once'));
end

reduced_files=all_files(keep_file);

end

function primary_variable=load_primary_variable_from_mat_file(file_name)

%Load first variable from MAT file
loaded_file=load(file_name);
loaded_variable_names=fieldnames(loaded_file);
if isempty(loaded_variable_names), error(['MAT file contains no variables: ' file_name]); end
primary_variable=loaded_file.(loaded_variable_names{1});

end

function value=get_structure_field_or_default(input_structure,field_name,default_value)

%Read structure field or default
if isstruct(input_structure) && isfield(input_structure,field_name) && ~isempty(input_structure.(field_name))
    value=input_structure.(field_name);
else
    value=default_value;
end

end

function value=get_scalar_or_nan(input_structure,field_name)

%Read scalar value
if isstruct(input_structure) && isfield(input_structure,field_name) && ~isempty(input_structure.(field_name))
    value=input_structure.(field_name);
    if numel(value)>1, value=value(1); end
else
    value=NaN;
end

end

function value=get_indexed_value_or_nan(input_structure,field_name,index_value)

%Read indexed value
if isstruct(input_structure) && isfield(input_structure,field_name) && numel(input_structure.(field_name))>=index_value
    value=input_structure.(field_name)(index_value);
else
    value=NaN;
end

end

function value=get_indexed_or_scalar_value_or_nan(input_structure,field_name,index_value)

%Read indexed or scalar value
if isstruct(input_structure) && isfield(input_structure,field_name) && ~isempty(input_structure.(field_name))
    field_value=input_structure.(field_name);
    if isscalar(field_value)
        value=field_value;
    elseif numel(field_value)>=index_value
        value=field_value(index_value);
    else
        value=NaN;
    end
else
    value=NaN;
end

end
