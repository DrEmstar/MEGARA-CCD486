function write_fits_output(objectname,final_data,barycentric_correction_flag,radial_velocity_measurement,calculate_lsd_profile)

%Import MATLAB FITS package
import matlab.io.*

%Store starting directory
starting_directory=pwd;
cleanup_directory=onCleanup(@() cd(starting_directory));

%Get target names
target_names=get_target_names_from_objectname(objectname);

%Write FITS output for each target
for target_index=1:numel(target_names)
    target_name=target_names{target_index};
    final_data_field_name=matlab.lang.makeValidName(target_name);
    if ~isfield(final_data,final_data_field_name) || isempty(final_data.(final_data_field_name)), continue; end
    target_directory=fullfile(starting_directory,target_name);
    if ~isfolder(target_directory), warning(['Reduced target directory does not exist: ' target_directory]); continue; end
    cd(target_directory)

    %Get target data
    target_data=final_data.(final_data_field_name);

    %Write regularised LSD profiles in FITS format
    if calculate_lsd_profile==1
        write_regularised_lsd_fits_file(target_name);
    end

    %decide if this is order-by-order
    target_data_field_names=fieldnames(target_data);
    if any(startsWith(target_data_field_names,'int'))
        %order-by-order data
        write_order_by_order_fits_files(target_name,target_data,barycentric_correction_flag);
    elseif any(strcmp(target_data_field_names,'fullwave')) && any(strcmp(target_data_field_names,'fullint'))
        %full merged spectra
        write_merged_fits_files(target_name,target_data,barycentric_correction_flag,radial_velocity_measurement);
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

function write_regularised_lsd_fits_file(target_name)

%Write regularised LSD profiles in FITS format
regularised_lsd_file_name=['regularised_lsd_profile_' target_name '.mat'];
if ~isfile(regularised_lsd_file_name), warning(['Regularised LSD profile file not found: ' regularised_lsd_file_name]); return; end

loaded_lsd_file=load(regularised_lsd_file_name);
if ~isfield(loaded_lsd_file,'regularised_lsd_profile'), warning(['Regularised LSD profile file does not contain regularised_lsd_profile: ' regularised_lsd_file_name]); return; end

regularised_lsd_profile=loaded_lsd_file.regularised_lsd_profile;
if ~isfield(regularised_lsd_profile,'velocity') || ~isfield(regularised_lsd_profile,'intensity'), warning(['regularised_lsd_profile is missing velocity or intensity in ' regularised_lsd_file_name]); return; end

[velocity,intensity]=standardise_lsd_profile_arrays(regularised_lsd_profile.velocity,regularised_lsd_profile.intensity,regularised_lsd_file_name);

%Output image columns are: velocity, profile 1, profile 2, ...
regularised_lsd_image=cat(2,velocity(:),intensity');

fits_lsd_filename=['regularised_lsd_profile_' target_name '.fits'];
file_pointer=create_or_replace_fits_file(fits_lsd_filename);

matlab.io.fits.createImg(file_pointer,'double_img',size(regularised_lsd_image));
matlab.io.fits.writeImg(file_pointer,regularised_lsd_image);

matlab.io.fits.writeKey(file_pointer,'PRODUCT','LSD','Line-profile product type');
matlab.io.fits.writeKey(file_pointer,'LSDMETH','REGINV','Regularised inverse-problem LSD');
matlab.io.fits.writeKey(file_pointer,'COL1','VELOCITY','First column is velocity in km/s');

if isfield(regularised_lsd_profile,'normalisation')
    write_optional_fits_key(file_pointer,'LSDNORM',char(string(regularised_lsd_profile.normalisation)),'LSD profile normalisation');
end

if isfield(regularised_lsd_profile,'velocity_step')
    write_optional_fits_key(file_pointer,'VELSTEP',regularised_lsd_profile.velocity_step,'LSD velocity step');
end

if isfield(regularised_lsd_profile,'regularisation_lambda')
    write_optional_fits_key(file_pointer,'REGLAM',regularised_lsd_profile.regularisation_lambda,'LSD regularisation lambda');
end

if isfield(regularised_lsd_profile,'regularisation_order')
    write_optional_fits_key(file_pointer,'REGORD',regularised_lsd_profile.regularisation_order,'LSD regularisation order');
end

if isfield(regularised_lsd_profile,'number_mask_lines')
    write_optional_fits_key(file_pointer,'NLINES',regularised_lsd_profile.number_mask_lines,'Number of LSD mask lines');
end

matlab.io.fits.closeFile(file_pointer);

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

function write_order_by_order_fits_files(target_name,target_data,barycentric_correction_flag)

%Find reduced files
file_list=list_reduced_spectrum_files();
if isempty(file_list), warning(['No reduced J*.mat files found for ' target_name]); return; end

%Find available orders
order_numbers=get_available_order_numbers(target_data);
if isempty(order_numbers), warning(['No order-by-order data found for ' target_name]); return; end
order_numbers=sort(order_numbers,'descend');

%Write one FITS file per observation
for observation_index=1:length(file_list)
    fits_file_root=file_list(observation_index).name(1:end-4);
    reduced_frame=load_primary_variable_from_mat_file(file_list(observation_index).name);
    fits_image=[];
    for order_index=1:numel(order_numbers)
        order_number=order_numbers(order_index);
        wave_field_name=['wave' num2str(order_number)];
        intensity_field_name=['int' num2str(order_number)];
        weights_field_name=['weights' num2str(order_number)];
        if ~isfield(target_data,wave_field_name) || ~isfield(target_data,intensity_field_name) || ~isfield(target_data,weights_field_name), continue; end
        if size(target_data.(intensity_field_name),1)<observation_index || size(target_data.(weights_field_name),1)<observation_index, continue; end
        wavelength_axis=target_data.(wave_field_name)(:);
        intensity_axis=target_data.(intensity_field_name)(observation_index,:)';
        weights_axis=target_data.(weights_field_name)(observation_index,:)';
        order_axis=order_number.*ones(length(wavelength_axis),1);
        fits_file_part=cat(2,wavelength_axis,intensity_axis,weights_axis,order_axis);
        fits_image=cat(1,fits_image,fits_file_part);
    end

    %header
    header_values=build_header_values(fits_file_root,reduced_frame,barycentric_correction_flag);

    %write to fits file
    fits_filename=[fits_file_root '_reduced.fits'];
    file_pointer=create_or_replace_fits_file(fits_filename);

    %image
    matlab.io.fits.createImg(file_pointer,'double_img',size(fits_image));
    matlab.io.fits.writeImg(file_pointer,fits_image);

    %header
    write_common_header_keys(file_pointer,header_values);
    matlab.io.fits.closeFile(file_pointer);
end

end

function write_merged_fits_files(target_name,target_data,barycentric_correction_flag,radial_velocity_measurement)

%Find reduced files
file_list=list_reduced_spectrum_files();
if isempty(file_list), warning(['No reduced J*.mat files found for ' target_name]); return; end

%Write one FITS file per observation
for observation_index=1:length(file_list)
    fits_file_root=file_list(observation_index).name(1:end-4);
    reduced_frame=load_primary_variable_from_mat_file(file_list(observation_index).name);
    fits_image=[];

    %wavelength axis
    fits_image(:,1)=target_data.fullwave(:);

    %intensity axis
    fits_image(:,2)=target_data.fullint(observation_index,:)';

    %header
    header_values=build_header_values(fits_file_root,reduced_frame,barycentric_correction_flag);
    if radial_velocity_measurement==1
        header_values.systemic_velocity=get_scalar_or_nan(target_data,'systemic_velocity');
        header_values.radial_velocity=get_indexed_value_or_nan(target_data,'radial_velocity',observation_index);
        header_values.radial_velocity_error=get_indexed_or_scalar_value_or_nan(target_data,'radial_velocity_error',observation_index);
        header_values.vbroad=get_broadening_value(target_data,'vbroad','vsini',observation_index);
        header_values.vbroad_error=get_broadening_value(target_data,'vbroad_error','vsini_error',observation_index);
        header_values.vsini=get_indexed_value_or_nan(target_data,'vsini',observation_index);
        header_values.vsini_error=get_indexed_or_scalar_value_or_nan(target_data,'vsini_error',observation_index);
    end

    %write to fits file
    fits_filename=[fits_file_root '_reduced.fits'];
    file_pointer=create_or_replace_fits_file(fits_filename);

    %image
    matlab.io.fits.createImg(file_pointer,'double_img',size(fits_image));
    matlab.io.fits.writeImg(file_pointer,fits_image);

    %header
    write_common_header_keys(file_pointer,header_values);
    if radial_velocity_measurement==1
        write_radial_velocity_header_keys(file_pointer,header_values);
    end
    matlab.io.fits.closeFile(file_pointer);
end

end

function header_values=build_header_values(fits_file_root,reduced_frame,barycentric_correction_flag)

%Initialise header values
header_values=struct();
header_values.filename=fits_file_root;
header_values.jd=get_structure_field_or_default(reduced_frame,'jd',NaN);

%Read header metadata if FITS file is available
fits_file_name=[fits_file_root '.fit'];
if isfile(fits_file_name)
    header=get_header(fits_file_name);
    headers=get_required_headers_from_header(header);
    header_values.date=get_structure_field_or_default(headers,'DATE_OBS',get_structure_field_or_default(reduced_frame,'expdate',''));
    header_values.exposure_start=get_structure_field_or_default(headers,'REC_STRT','');
    header_values.exposure_time=get_structure_field_or_default(headers,'EXPTIME',get_structure_field_or_default(reduced_frame,'exptime',NaN));
    header_values.flux_weighted_midtime=get_structure_field_or_default(headers,'HERCFWMT',NaN);
else
    header_values.date=get_structure_field_or_default(reduced_frame,'expdate','');
    header_values.exposure_start=get_structure_field_or_default(reduced_frame,'midtime_utc','');
    header_values.exposure_time=get_structure_field_or_default(reduced_frame,'exptime',NaN);
    header_values.flux_weighted_midtime=NaN;
end

%Read stored reduction metadata
header_values.barycentric_correction=get_structure_field_or_default(reduced_frame,'bcorr',NaN);
header_values.barycentric_correction_applied=barycentric_correction_flag;
header_values.signal_to_noise=get_structure_field_or_default(reduced_frame,'signal_to_noise',NaN);

end

function write_common_header_keys(file_pointer,header_values)

%Write common FITS header keys
matlab.io.fits.writeKey(file_pointer,'FILENAME',header_values.filename,'Original filename');
write_optional_fits_key(file_pointer,'JD',header_values.jd,'Reduced heliocentric mid-obs JD');
write_optional_fits_key(file_pointer,'DATE',header_values.date,'Observation date');
write_optional_fits_key(file_pointer,'START',header_values.exposure_start,'Exposure start time UTC');
write_optional_fits_key(file_pointer,'EXPTIME',header_values.exposure_time,'Exposure duration');
write_optional_fits_key(file_pointer,'FWMT',header_values.flux_weighted_midtime,'Flux-weighted exposure midtime UTC');
write_optional_fits_key(file_pointer,'BCORR',header_values.barycentric_correction,'Barycentric correction velocity');
write_optional_fits_key(file_pointer,'BCAPP',header_values.barycentric_correction_applied,'Barycentric correction applied');
write_optional_fits_key(file_pointer,'SNR',header_values.signal_to_noise,'Signal-to-noise');

end

function write_radial_velocity_header_keys(file_pointer,header_values)

%Write radial-velocity FITS header keys
write_optional_fits_key(file_pointer,'SYSVEL',header_values.systemic_velocity,'Systemic velocity');
write_optional_fits_key(file_pointer,'RV',header_values.radial_velocity,'Radial velocity');
write_optional_fits_key(file_pointer,'RVERR',header_values.radial_velocity_error,'Radial velocity error');
write_optional_fits_key(file_pointer,'VSINI',header_values.vsini, ...
    'Instrument-corrected width; macroturbulence retained');
write_optional_fits_key(file_pointer,'VSIERR',header_values.vsini_error,'Vsini error');
write_optional_fits_key(file_pointer,'VBROAD',header_values.vbroad, ...
    'Observed width including instrumental broadening');
write_optional_fits_key(file_pointer,'VBRDERR',header_values.vbroad_error,'Vbroad error');

end

function value=get_broadening_value(target_data,preferred_field,legacy_field,observation_index)

%Read the new field, falling back to the legacy alias for older products.
if isfield(target_data,preferred_field)
    value=get_indexed_or_scalar_value_or_nan(target_data,preferred_field,observation_index);
else
    value=get_indexed_or_scalar_value_or_nan(target_data,legacy_field,observation_index);
end

end

function write_optional_fits_key(file_pointer,key_name,value,comment)

%CFITSIO cannot write NaN or Inf numeric keyword values. Missing optional
%metadata are represented internally by NaN or empty values, so omit those
%keywords while retaining every valid value.
if isempty(value), return; end
if isnumeric(value) || islogical(value)
    if ~isscalar(value) || ~isfinite(double(value)), return; end
elseif isstring(value)
    if ~isscalar(value) || strlength(value)==0, return; end
    value=char(value);
elseif ~ischar(value)
    return
end
matlab.io.fits.writeKey(file_pointer,key_name,value,comment);

end

function file_pointer=create_or_replace_fits_file(fits_filename)

%Create or replace FITS file
if isfile(fits_filename), delete(fits_filename); end
file_pointer=matlab.io.fits.createFile(fits_filename);

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
