%BACKFILL_BCORR Add missing barycentric corrections to reduced MEGARA frames.
% Edit objectname below, then run this script with MEGARA 2.1 on the MATLAB
% path. Files are updated in <target>/reduced frames using the MATLAB-only
% barycentric correction. Existing finite corrections are left unchanged.

clearvars
objectname='hd_39060';

reduced_data_pointer=which('pointer_reduced_data_directory.m');
if isempty(reduced_data_pointer)
    error('pointer_reduced_data_directory.m is not on the MATLAB path.')
end
reduced_data_directory=fileparts(reduced_data_pointer);
frames_directory=fullfile(reduced_data_directory,objectname,'reduced frames');
final_files=dir(fullfile(frames_directory,'J*.mat'));
if isempty(final_files)
    error('No reduced frames were found in %s.',frames_directory)
end

for file_index=1:numel(final_files)
    frame_path=fullfile(final_files(file_index).folder,final_files(file_index).name);
    loaded=load(frame_path);
    variable_names=fieldnames(loaded);
    if numel(variable_names)~=1 || ~isstruct(loaded.(variable_names{1}))
        warning('Skipping %s: expected one reduced-frame structure.',frame_path)
        continue
    end
    variable_name=variable_names{1};
    frame=loaded.(variable_name);
    if isfield(frame,'bcorr') && isscalar(frame.bcorr) && isfinite(frame.bcorr)
        continue
    end
    if ~isfield(frame,'midtime_tt')
        midtime_utc=datetime([strip_quotes(frame.expdate) ' ' ...
            strip_quotes(frame.midtime_hrsp)], ...
            'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        frame.midtime_utc=['"' char(midtime_utc,'HH:mm:ss.SSS') '"'];
        frame.midtime_tt=['"' char(midtime_utc+seconds(69.184), ...
            'HH:mm:ss.SSS') '"'];
    end
    fits_name=[variable_name '.fit'];
    correction=get_barycentric_correction( ...
        frame.starname,frame.expdate,frame.midtime_tt,fits_name);
    frame.bcorr=str2double(correction);
    if ~isfinite(frame.bcorr)
        warning('No finite barycentric correction was returned for %s.',frame_path)
        continue
    end
    output=struct(variable_name,frame);
    save(frame_path,'-struct','output')
end

function value=strip_quotes(value)
value=strrep(char(value),'"','');
end
