function [final_data,objectname,star_database]=megara_control_file(options)
warning on
warning('off','backtrace')

%Ensure the bundled MATLAB barycentric correction is available even when
%this subfolder was added after the user's saved MATLAB path was created.
reduction_code_directory=fileparts(mfilename('fullpath'));
matlab_barycentric_directory=fullfile(reduction_code_directory, ...
    'barycentric_correction','matlab_barycentric');
if ~isfolder(matlab_barycentric_directory)
    error('Bundled MATLAB barycentric-correction folder not found: %s', ...
        matlab_barycentric_directory)
end
addpath(matlab_barycentric_directory,'-begin')

options=validate_megara_options(options);

final_data=struct();
star_database=[];
found_objectnames=struct();

home_directory_megara_file=which('pointer_file_megara.mat');
if isempty(home_directory_megara_file), error('Could not find pointer_file_megara.mat on the MATLAB path.'); end
home_directory_megara=fileparts(home_directory_megara_file);
original_directory=pwd;
cleanup_object=onCleanup(@() cd(original_directory));
cd(home_directory_megara)

[raw_data_location,reduced_data_location]=define_directory_locations;

objectname=options.objectname;

if ischar(objectname) && ~isempty(objectname) && ~any(strcmpi(objectname,{'musician','tess','exocomet'}))
    single_name=objectname;
    objectname=struct();
    objectname.target1=single_name;
end

if options.star_census
    found_objectnames=create_observation_summaries(raw_data_location,reduced_data_location,options.reduction_folders,objectname);
end

if isstruct(objectname)

    target_names=get_target_names_from_object_structure(objectname);

    for target_index=1:numel(target_names)

        single_objectname=canonicalise_single_objectname(target_names{target_index});
        disp(['Reducing object ' single_objectname])

        star_database=load_star_database_for_object(single_objectname);

        if isempty(star_database) && options.star_census
            star_database=load_star_database_for_object(single_objectname);
        end

        if isempty(star_database)
            star_database.exptime=[];
        end

        [single_object_teff,single_object_logg,single_object_vsini]=get_stellar_parameters_for_object(single_objectname,options);
        final_data_single_object=reduce_single_object(single_objectname,single_object_teff,single_object_logg,single_object_vsini,options,raw_data_location,reduced_data_location);
        drawnow
        final_data=store_final_data(final_data,single_objectname,final_data_single_object);
        save_outputs_for_object(single_objectname,final_data_single_object,options,reduced_data_location);

    end

elseif ischar(objectname) && any(strcmpi(objectname,{'musician','tess','exocomet'}))

    database_mode=objectname;
    database=load_target_database(database_mode);
    output_objectname=struct();
    number_of_objects=0;

    for database_index=setdiff(1:size(database,1),[46,48,59])

        single_objectname=canonicalise_single_objectname(database{database_index,1});
        single_object_teff=database{database_index,2};
        single_object_logg=database{database_index,3};
        single_object_vsini=database{database_index,4};

        number_of_objects=number_of_objects+1;
        output_objectname.(['target' num2str(number_of_objects)])=single_objectname;

        disp(['Reducing object ' single_objectname])

        star_database=load_star_database_for_object(single_objectname);
        if isempty(star_database)
            star_database.exptime=[];
        end

        final_data_single_object=reduce_single_object(single_objectname,single_object_teff,single_object_logg,single_object_vsini,options,raw_data_location,reduced_data_location);
        drawnow
        final_data=store_final_data(final_data,single_objectname,final_data_single_object);
        save_outputs_for_object(single_objectname,final_data_single_object,options,reduced_data_location);

    end

    objectname=output_objectname;

elseif isempty(objectname)

    disp('Reducing all objects')

    if ~options.star_census
        error('options.objectname=[] requires options.star_census=true, because the code needs a census to discover all target names.')
    end

    if isempty(fieldnames(found_objectnames))
        error('No objects were found by create_observation_summaries. Check raw data folders and options.reduction_folders.')
    end

    objectname=canonicalise_objectname_input(found_objectnames);
    target_names=get_target_names_from_object_structure(objectname);

    fprintf('A total of %.0f objects will be reduced\n',numel(target_names))

    for target_index=1:numel(target_names)

        single_objectname=canonicalise_single_objectname(target_names{target_index});
        disp(['Reducing object ' single_objectname])

        [single_object_teff,single_object_logg,single_object_vsini]=get_stellar_parameters_for_object(single_objectname,options);
        final_data_single_object=reduce_single_object(single_objectname,single_object_teff,single_object_logg,single_object_vsini,options,raw_data_location,reduced_data_location);
        drawnow
        final_data=store_final_data(final_data,single_objectname,final_data_single_object);
        save_outputs_for_object(single_objectname,final_data_single_object,options,reduced_data_location);

    end

    star_database='See individual database files';

else

    error('options.objectname must be empty, a target structure, a single target name, or one of musician, tess, or exocomet.');

end

temporary_files=dir('temp*');
if ~isempty(temporary_files), delete('temp*'); end

cd(reduced_data_location)

if options.ascii_output
    write_ascii_output(objectname,final_data,options.apply_barycentric_correction,options.radial_velocity_measurement,options.calculate_lsd_profile)
end

if options.fits_output
    write_fits_output(objectname,final_data,options.apply_barycentric_correction,options.radial_velocity_measurement,options.calculate_lsd_profile)
end

%A failed run returns to the directory from which MEGARA was started via
%cleanup_object.  After a successful run, deliberately finish beside the
%products: in the target directory for one target, or at the reduced-data
%root when more than one target was requested.
completed_run_directory=get_completed_run_directory(objectname,reduced_data_location);
delete(cleanup_object)
cd(completed_run_directory)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function completed_run_directory=get_completed_run_directory(objectname,reduced_data_location)

completed_run_directory=reduced_data_location;

if ~isstruct(objectname)
    return
end

target_names=get_target_names_from_object_structure(objectname);
if numel(target_names)~=1
    return
end

single_target_directory=fullfile(reduced_data_location, ...
    canonicalise_single_objectname(target_names{1}));
if isfolder(single_target_directory)
    completed_run_directory=single_target_directory;
else
    warning(['The single-target output directory was not found. ' ...
        'MEGARA will finish in the reduced-data root: ' ...
        single_target_directory])
end

end

function options=validate_megara_options(options)

if ~isstruct(options), error('Input must be an options structure.'); end

if isfield(options,['cross' '_correlation'])
    error('Obsolete line-profile option found. Use options.calculate_lsd_profile only.')
end

if isfield(options,'max_stellar_spectra') || isfield(options,'stellar_spectrum_start_index')
    error('options.max_stellar_spectra and options.stellar_spectrum_start_index are obsolete for production post-reduction. Use options.chunked_post_reduction and options.post_reduction_chunk_size instead.')
end

if ~isfield(options,'chunked_post_reduction') || isempty(options.chunked_post_reduction)
    options.chunked_post_reduction=false;
end

if ~isfield(options,'post_reduction_chunk_size') || isempty(options.post_reduction_chunk_size)
    options.post_reduction_chunk_size=500;
end

if ~isfield(options,'extended_wavelength_range') || ...
        isempty(options.extended_wavelength_range)
    options.extended_wavelength_range=false;
end

if ~isfield(options,'skip_manual_file_check') || ...
        isempty(options.skip_manual_file_check)
    options.skip_manual_file_check=false;
end

required_fields={'objectname','reduction_folders','teff','logg','vsini','star_census','skip_reduction','skip_post_reduction','blue_data_chop_value','apply_barycentric_correction','apply_continuum_fit','manual_continuum_fit','template_continuum_fit','automatic_continuum_fit','order_merge','overwrite_full_data','overwrite_post_reduction','calculate_lsd_profile','radial_velocity_measurement','extend_velocities','suppress_figures','ascii_output','fits_output','FAMIAS_output'};

for required_field_index=1:numel(required_fields)
    if ~isfield(options,required_fields{required_field_index}), error(['Missing options field: options.' required_fields{required_field_index}]); end
end

if isstring(options.objectname), options.objectname=char(options.objectname); end
if ~(ischar(options.objectname) || isstruct(options.objectname) || isempty(options.objectname)), error('options.objectname must be a character vector, string scalar, structure, or empty.'); end

options.objectname=canonicalise_objectname_input(options.objectname);

if ischar(options.reduction_folders), options.reduction_folders={options.reduction_folders}; end
if isstring(options.reduction_folders), options.reduction_folders=cellstr(options.reduction_folders); end
if ~(iscellstr(options.reduction_folders) || isempty(options.reduction_folders)), error('options.reduction_folders must be a cell array of character vectors, a string array, a character vector, or empty.'); end

if ~isempty(options.teff) && ~(isnumeric(options.teff) && isscalar(options.teff) && isfinite(options.teff)), error('options.teff must be a finite scalar numeric value or empty.'); end
if ~isempty(options.logg) && ~(isnumeric(options.logg) && isscalar(options.logg) && isfinite(options.logg)), error('options.logg must be a finite scalar numeric value or empty.'); end
if ~isempty(options.vsini) && ~(isnumeric(options.vsini) && isscalar(options.vsini)), error('options.vsini must be a scalar numeric value or empty.'); end

if ~(isnumeric(options.blue_data_chop_value) && isscalar(options.blue_data_chop_value) && options.blue_data_chop_value>=0 && isfinite(options.blue_data_chop_value)), error('options.blue_data_chop_value must be a finite scalar number >= 0.'); end

if ~(isnumeric(options.post_reduction_chunk_size) && isscalar(options.post_reduction_chunk_size) && isfinite(options.post_reduction_chunk_size) && options.post_reduction_chunk_size>=1 && floor(options.post_reduction_chunk_size)==options.post_reduction_chunk_size)
    error('options.post_reduction_chunk_size must be a positive integer.')
end

logical_fields={'star_census','skip_reduction','skip_post_reduction','skip_manual_file_check','chunked_post_reduction','extended_wavelength_range','apply_barycentric_correction','apply_continuum_fit','manual_continuum_fit','template_continuum_fit','automatic_continuum_fit','order_merge','overwrite_full_data','overwrite_post_reduction','calculate_lsd_profile','radial_velocity_measurement','extend_velocities','suppress_figures','ascii_output','fits_output','FAMIAS_output'};

for logical_field_index=1:numel(logical_fields)
    logical_field_name=logical_fields{logical_field_index};
    if ~(islogical(options.(logical_field_name)) && isscalar(options.(logical_field_name))) && ~(isnumeric(options.(logical_field_name)) && isscalar(options.(logical_field_name)) && any(options.(logical_field_name)==[0 1])), error(['options.' logical_field_name ' must be true or false.']); end
    options.(logical_field_name)=logical(options.(logical_field_name));
end

if options.order_merge && ~options.apply_continuum_fit, error('options.order_merge=true requires options.apply_continuum_fit=true.'); end
if sum([options.manual_continuum_fit options.template_continuum_fit ...
        options.automatic_continuum_fit])>1
    error('Select only one manual, template, or automatic continuum-fitting option.')
end
if options.automatic_continuum_fit && ~options.order_merge
    error('options.automatic_continuum_fit=true requires options.order_merge=true for the merged-spectrum SUPPNet stage.')
end
if options.template_continuum_fit && ~options.order_merge
    error('options.template_continuum_fit=true requires options.order_merge=true.')
end
if options.calculate_lsd_profile && ~options.order_merge, error('options.calculate_lsd_profile=true requires options.order_merge=true.'); end
if options.radial_velocity_measurement && ~options.calculate_lsd_profile, error('options.radial_velocity_measurement=true requires options.calculate_lsd_profile=true.'); end
if options.FAMIAS_output && (~options.overwrite_post_reduction || ~options.calculate_lsd_profile), error('options.FAMIAS_output=true requires options.overwrite_post_reduction=true and options.calculate_lsd_profile=true.'); end
if options.skip_post_reduction && (options.calculate_lsd_profile || options.radial_velocity_measurement || options.FAMIAS_output), warning('options.skip_post_reduction=true means regularised LSD profiles, radial velocity measurement, and FAMIAS output will probably not be regenerated.'); end

end

function objectname=canonicalise_objectname_input(objectname)

if isempty(objectname)
    return
end

if ischar(objectname)
    if any(strcmpi(objectname,{'musician','tess','exocomet'}))
        objectname=lower(objectname);
    else
        objectname=lower(objectname);
    end
    return
end

if isstruct(objectname)
    target_field_names=fieldnames(objectname);
    for target_index=1:numel(target_field_names)
        current_value=objectname.(target_field_names{target_index});
        if isstring(current_value), current_value=char(current_value); end
        if ~ischar(current_value), error('Each target name in options.objectname must be a character vector or string scalar.'); end
        objectname.(target_field_names{target_index})=canonicalise_single_objectname(current_value);
    end
end

end

function single_objectname=canonicalise_single_objectname(single_objectname)

if isstring(single_objectname)
    single_objectname=char(single_objectname);
end

if ~ischar(single_objectname)
    error('Object name must be a character vector or string scalar.')
end

single_objectname=lower(single_objectname);

end

function target_names=get_target_names_from_object_structure(objectname)

target_field_names=fieldnames(objectname);
number_of_targets=numel(target_field_names);
target_names=cell(number_of_targets,1);

for target_index=1:number_of_targets

    expected_target_field_name=['target' num2str(target_index)];

    if isfield(objectname,expected_target_field_name)
        target_names{target_index}=objectname.(expected_target_field_name);
    else
        target_names{target_index}=objectname.(target_field_names{target_index});
        warning(['Expected objectname.' expected_target_field_name ' but found objectname.' target_field_names{target_index} '. Using available field order instead.']);
    end

    if isstring(target_names{target_index}), target_names{target_index}=char(target_names{target_index}); end
    if ~ischar(target_names{target_index}), error('Each target name in options.objectname must be a character vector or string scalar.'); end
    target_names{target_index}=canonicalise_single_objectname(target_names{target_index});

end

end

function [teff,logg,vsini]=get_stellar_parameters_for_object(single_objectname,options)

teff=options.teff;
logg=options.logg;
vsini=options.vsini;

if is_finite_numeric_scalar(teff) && is_finite_numeric_scalar(logg) && is_finite_numeric_scalar(vsini)
    return
end

database_names={'musician','tess','exocomet'};

for database_index=1:numel(database_names)

    database=load_target_database(database_names{database_index});
    object_index=find(strcmpi(database(:,1),single_objectname),1,'first');

    if isempty(object_index)
        continue
    end

    if ~is_finite_numeric_scalar(teff)
        teff=database{object_index,2};
    end

    if ~is_finite_numeric_scalar(logg)
        logg=database{object_index,3};
    end

    if ~is_finite_numeric_scalar(vsini)
        vsini=database{object_index,4};
    end

    break

end

if ~is_finite_numeric_scalar(teff)
    error(['teff is missing for ' single_objectname '. Set options.teff or add the target to the MUSICIAN/TESS/exocomet database.'])
end

if ~is_finite_numeric_scalar(logg)
    error(['logg is missing for ' single_objectname '. Set options.logg or add the target to the MUSICIAN/TESS/exocomet database.'])
end

if ~is_finite_numeric_scalar(vsini)
    vsini=NaN;
end

end

function flag=is_finite_numeric_scalar(value)

flag=isnumeric(value) && isscalar(value) && isfinite(value);

end

function database=load_target_database(objectname)

if strcmpi(objectname,'musician')
    loaded_database=load('mus_database_cool.mat','mus_database');
    database=loaded_database.mus_database;
elseif strcmpi(objectname,'tess')
    loaded_database=load('tess_database.mat','tess_database');
    database=loaded_database.tess_database;
elseif strcmpi(objectname,'exocomet')
    loaded_database=load('exocomet_database.mat','exocomet_database');
    database=loaded_database.exocomet_database;
else
    error(['Unknown target database: ' objectname]);
end

end

function star_database=load_star_database_for_object(single_objectname)

star_database=[];
primary_database_file_name=['star_database_' single_objectname '.mat'];
primary_database_variable_name=['star_database_' single_objectname];

if isfile(primary_database_file_name)
    loaded_database=load(primary_database_file_name);
    star_database=get_database_variable_from_loaded_file(loaded_database,primary_database_variable_name,primary_database_file_name);
    return
end

if numel(single_objectname)>2

    secondary_database_file_name=['star_database_HD' single_objectname(3:end) '.mat'];
    secondary_database_variable_name=['star_database_HD' single_objectname(3:end)];

    if isfile(secondary_database_file_name)
        loaded_database=load(secondary_database_file_name);
        star_database=get_database_variable_from_loaded_file(loaded_database,secondary_database_variable_name,secondary_database_file_name);
    end

end

end

function star_database=get_database_variable_from_loaded_file(loaded_database,expected_variable_name,database_file_name)

if isfield(loaded_database,expected_variable_name)
    star_database=loaded_database.(expected_variable_name);
elseif isfield(loaded_database,'star_database')
    star_database=loaded_database.star_database;
else
    loaded_variable_names=fieldnames(loaded_database);
    star_database=loaded_database.(loaded_variable_names{1});
    warning(['Could not find expected variable ' expected_variable_name ' in ' database_file_name '. Using first variable in file instead.']);
end

end

function final_data_single_object=reduce_single_object(single_objectname,teff,logg,vsini,options,raw_data_location,reduced_data_location)

final_data_single_object=reduce_process_1_target(single_objectname,options,raw_data_location,reduced_data_location,teff,logg,vsini);

end

function final_data=store_final_data(final_data,single_objectname,final_data_single_object)

single_object_field_name=matlab.lang.makeValidName(single_objectname);

if ~strcmp(single_object_field_name,single_objectname)
    warning(['Object name ' single_objectname ' is not a valid MATLAB structure field name. Stored as final_data.' single_object_field_name '.']);
end

final_data.(single_object_field_name)=final_data_single_object;

end

function save_outputs_for_object(single_objectname,final_data_single_object,options,reduced_data_location)

if options.skip_post_reduction || isempty(final_data_single_object), return; end

reduced_target_directory=fullfile(reduced_data_location,single_objectname);
if ~isfolder(reduced_target_directory), error(['Reduced target directory does not exist: ' reduced_target_directory]); end

finaldata=final_data_single_object;

if options.order_merge
    final_data_file_name=fullfile(reduced_target_directory,['final_data_' single_objectname '.mat']);
else
    final_data_file_name=fullfile(reduced_target_directory,['final_data_unmerged_' single_objectname '.mat']);
end

save(final_data_file_name,'finaldata','-v7.3')

if isfield(finaldata,'regularised_lsd_profile')
    regularised_lsd_profile=finaldata.regularised_lsd_profile;
    regularised_lsd_profile_file=fullfile(reduced_target_directory,['regularised_lsd_profile_' single_objectname '.mat']);
    save(regularised_lsd_profile_file,'regularised_lsd_profile','-v7.3')
end

if options.FAMIAS_output && isfield(finaldata,'regularised_lsd_profile')
    current_directory=pwd;
    cleanup_directory=onCleanup(@() cd(current_directory));
    cd(reduced_target_directory)
    write_regularised_lsd_FAMIAS_output(single_objectname,finaldata)
end

end

function write_regularised_lsd_FAMIAS_output(objectname,finaldata)

if ~isfield(finaldata,'fullint'), error('FAMIAS output requires finaldata.fullint.'); end
if ~isfield(finaldata,'jd'), error('FAMIAS output requires finaldata.jd.'); end
if ~isfield(finaldata,'regularised_lsd_profile'), error('FAMIAS output requires finaldata.regularised_lsd_profile.'); end

regularised_lsd_profile=finaldata.regularised_lsd_profile;
number_noise_pixels=min(10000,size(finaldata.fullint,2));
intfilt=medfilt1(finaldata.fullint(:,1:number_noise_pixels)',11)';
resid=finaldata.fullint(:,1:number_noise_pixels)-intfilt;
noise_sigma=std(resid,[],2,'omitnan');
valid_noise=isfinite(noise_sigma) & noise_sigma>0;
if ~any(valid_noise)
    error('Could not calculate finite FAMIAS weights from finaldata.fullint.')
end
noise_sigma(~valid_noise)=median(noise_sigma(valid_noise));
weight=1./noise_sigma;

if ~isfolder('Regularised_LSD_FAMIAS_files')
    mkdir('Regularised_LSD_FAMIAS_files')
end

current_directory=pwd;
cleanup_directory=onCleanup(@() cd(current_directory));
cd('Regularised_LSD_FAMIAS_files')
write_files_for_FAMIAS(finaldata.jd,regularised_lsd_profile.velocity,regularised_lsd_profile.intensity,weight,['times_' get_object_display_name(objectname) '.txt'],[get_object_display_name(objectname) '_'])

end

function object_display_name=get_object_display_name(objectname)

if numel(objectname)>=4
    object_display_name=objectname(4:end);
else
    object_display_name=objectname;
end

end
