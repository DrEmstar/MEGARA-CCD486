function finaldata=process_1_target(objectname,options,reduced_data_location,teff,logg,vsini)

%PROCESS_1_TARGET Run post-reduction for one target.

%Initialise output
finaldata=[];

%Normalise object name
if isstring(objectname)
    objectname=char(objectname);
end

if ~ischar(objectname) || isempty(objectname)
    error('objectname must be a non-empty character vector or string scalar.')
end

%Validate options needed by this stage
options=validate_process_options(options);

%Validate reduced-data location
if isstring(reduced_data_location)
    reduced_data_location=char(reduced_data_location);
end

if ~ischar(reduced_data_location) || isempty(reduced_data_location)
    error('reduced_data_location must be a non-empty character vector or string scalar.')
end

if ~isfolder(reduced_data_location)
    error(['Reduced data location does not exist: ' reduced_data_location])
end

%Move to exact reduced target folder before post-reduction
reduced_target_directory=fullfile(reduced_data_location,objectname);

if ~isfolder(reduced_target_directory)
    error(['Reduced target directory does not exist: ' reduced_target_directory])
end

starting_directory=pwd;
cleanup_directory=onCleanup(@() cd(starting_directory));
cd(reduced_target_directory)

%Report post-reduction directory
fprintf('Post-reduction target directory: %s\n',pwd)

%Check stellar parameters needed for regularised LSD line-mask generation
validate_stellar_parameters_for_lsd(options,teff,logg,vsini,objectname);

%Run post-reduction. post_reduction can return myflag=0 if it removed bad
%files and needs to be rerun on the cleaned file set.
maximum_post_reduction_attempts=2;
myflag=0;

for post_reduction_attempt=1:maximum_post_reduction_attempts

    fprintf('Post-reduction pass %.0f for %s\n',post_reduction_attempt,objectname)

    [myflag,finaldata]=post_reduction(objectname,options,teff,logg,vsini);

    if myflag==1
        break
    end

    if post_reduction_attempt==1
        disp('post_reduction removed bad files. Reprocessing once using the cleaned file set.')
    end

    drawnow

end

%Stop if cleanup did not converge after one reprocessing pass
if myflag~=1
    error(['post_reduction did not complete after one cleanup and reprocessing pass for ' objectname '. Bad-file cleanup is still removing files after the cleaned rerun. The cleanup criteria in post_reduction.m need to be made more complete before expensive processing starts.'])
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function options=validate_process_options(options)

%Check options structure
if ~isstruct(options)
    error('options must be a structure.')
end

%Check fields used directly by this function
required_fields={'calculate_lsd_profile','radial_velocity_measurement','FAMIAS_output','chunked_post_reduction','post_reduction_chunk_size'};

for required_field_index=1:numel(required_fields)
    if ~isfield(options,required_fields{required_field_index})
        error(['Missing options field: options.' required_fields{required_field_index}])
    end
end

%Reject removed MEGARA 1 compatibility option
if isfield(options,'cross_correlation')
    error('options.cross_correlation is obsolete in this no-backward-compatibility MEGARA 2.1 workflow. Use options.calculate_lsd_profile only.')
end

%Reject obsolete partial-processing options
if isfield(options,'max_stellar_spectra') || isfield(options,'stellar_spectrum_start_index') || isfield(options,'processing_observation_numbers')
    error('Partial post-reduction selection options are obsolete. Use options.chunked_post_reduction and options.post_reduction_chunk_size to process all spectra in memory-safe chunks.')
end

%Validate chunked post-reduction options
if ~(islogical(options.chunked_post_reduction) && isscalar(options.chunked_post_reduction)) && ~(isnumeric(options.chunked_post_reduction) && isscalar(options.chunked_post_reduction) && any(options.chunked_post_reduction==[0 1]))
    error('options.chunked_post_reduction must be true or false.')
end

options.chunked_post_reduction=logical(options.chunked_post_reduction);

if ~(isnumeric(options.post_reduction_chunk_size) && isscalar(options.post_reduction_chunk_size) && isfinite(options.post_reduction_chunk_size) && options.post_reduction_chunk_size>=1 && floor(options.post_reduction_chunk_size)==options.post_reduction_chunk_size)
    error('options.post_reduction_chunk_size must be a positive integer.')
end

end

function validate_stellar_parameters_for_lsd(options,teff,logg,vsini,objectname)

%Only regularised LSD generation needs synthetic-spectrum stellar parameters
needs_lsd_line_mask=options.calculate_lsd_profile==1 || options.radial_velocity_measurement==1 || options.FAMIAS_output==1;

if ~needs_lsd_line_mask
    return
end

%teff and logg must be finite because they determine the synthetic spectrum
if ~is_finite_numeric_scalar(teff)
    error(['teff must be a finite numeric scalar for regularised LSD line-mask generation for ' objectname '.'])
end

if ~is_finite_numeric_scalar(logg)
    error(['logg must be a finite numeric scalar for regularised LSD line-mask generation for ' objectname '.'])
end

%vsini can be missing because post_reduction uses a safe finite value for
%synthetic-spectrum line-mask generation only
if ~is_finite_numeric_scalar(vsini)
    warning(['vsini is missing or non-finite for ' objectname '. post_reduction should use vsini=0 km/s for synthetic-spectrum and LSD line-mask generation only.'])
end

end

function flag=is_finite_numeric_scalar(value)

%Check finite numeric scalar
flag=isnumeric(value) && isscalar(value) && isfinite(value);

end
