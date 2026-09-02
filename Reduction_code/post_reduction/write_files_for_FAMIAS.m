function write_files_for_FAMIAS(jd,profile_axis,profile_intensity,weights,times_name,filestubname)

%Normalise output names
if isstring(times_name), times_name=char(times_name); end
if isstring(filestubname), filestubname=char(filestubname); end

%Validate inputs
validate_FAMIAS_inputs(jd,profile_axis,profile_intensity,weights,times_name,filestubname);

%Force standard orientations
jd=jd(:);
profile_axis=profile_axis(:).';
weights=weights(:);

if size(profile_intensity,2)==numel(profile_axis)
    %Already spectra x profile axis
elseif size(profile_intensity,1)==numel(profile_axis)
    profile_intensity=profile_intensity.';
else
    error('profile_intensity has incompatible dimensions after validation.')
end

%Open times file
file_identifier=fopen(times_name,'w+');
if file_identifier==-1, error(['Could not open FAMIAS times file: ' times_name]); end
cleanup_object=onCleanup(@() fclose(file_identifier));

%Write profiles and times file
for spectrum_index=1:numel(jd)
    output_file_name=[filestubname num2str(spectrum_index) '.dat'];
    fprintf(file_identifier,'%s %.6f %.6f\n',output_file_name,jd(spectrum_index),weights(spectrum_index));
    tempvar=cat(2,profile_axis(:),profile_intensity(spectrum_index,:)');
    save(output_file_name,'tempvar','-ascii')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_FAMIAS_inputs(jd,profile_axis,profile_intensity,weights,times_name,filestubname)

%Check time values
if ~isnumeric(jd) || ~isvector(jd), error('jd must be a numeric vector.'); end
if any(~isfinite(jd(:))), error('jd must contain only finite values.'); end

%Check profile axis
if ~isnumeric(profile_axis) || ~isvector(profile_axis), error('profile_axis must be a numeric vector.'); end
if any(~isfinite(profile_axis(:))), error('profile_axis must contain only finite values.'); end
if numel(profile_axis)<2, error('profile_axis must contain at least two values.'); end

%Check profile intensity matrix
if ~isnumeric(profile_intensity) || ndims(profile_intensity)~=2, error('profile_intensity must be a numeric 2D matrix.'); end

if size(profile_intensity,1)==numel(jd) && size(profile_intensity,2)==numel(profile_axis)
    %Already spectra x profile axis
elseif size(profile_intensity,2)==numel(jd) && size(profile_intensity,1)==numel(profile_axis)
    %Profile axis x spectra; this is allowed and transposed in the main function
else
    error('profile_intensity must have one dimension matching jd and the other matching profile_axis.')
end

if any(~isfinite(profile_intensity(:))), error('profile_intensity must contain only finite values.'); end

%Check weights
if ~isnumeric(weights) || ~isvector(weights), error('weights must be a numeric vector.'); end
if numel(weights)~=numel(jd), error('weights must have the same number of elements as jd.'); end
if any(~isfinite(weights(:))), error('weights must contain only finite values.'); end

%Check output names
if ~ischar(times_name) || isempty(times_name), error('times_name must be a non-empty character vector or string scalar.'); end
if ~ischar(filestubname) || isempty(filestubname), error('filestubname must be a non-empty character vector or string scalar.'); end

end