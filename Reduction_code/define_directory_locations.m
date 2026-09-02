function [raw_data_location,reduced_data_location]=define_directory_locations()
%setting up where all the directories are

%Ignore the example pointer files supplied inside the active MEGARA tree.
megara_pointer_file=which('pointer_file_megara.mat');
if isempty(megara_pointer_file)
    error('Could not find pointer_file_megara.mat on the MATLAB path.')
end
megara_directory=fileparts(megara_pointer_file);

%Find raw data pointer file
raw_data_pointer_file=find_external_pointer_file( ...
    'pointer_raw_data_directory.m',megara_directory);
if isempty(raw_data_pointer_file)
    error(['No external raw data directory found. Copy ' ...
        'pointer_raw_data_directory.m to the raw-data root and add that ' ...
        'directory to the MATLAB path. Pointer copies inside the MEGARA ' ...
        'program directory are ignored.'])
end

%Set raw data directory
raw_data_location=fileparts(raw_data_pointer_file);

%Find reduced data pointer file
reduced_data_pointer_file=find_external_pointer_file( ...
    'pointer_reduced_data_directory.m',megara_directory);
if isempty(reduced_data_pointer_file)
    error(['No external reduced data directory found. Copy ' ...
        'pointer_reduced_data_directory.m to the reduced-data root and ' ...
        'add that directory to the MATLAB path. Pointer copies inside ' ...
        'the MEGARA program directory are ignored.'])
end

%Set reduced data directory
reduced_data_location=fileparts(reduced_data_pointer_file);

end

function pointer_file=find_external_pointer_file(pointer_name,megara_directory)

pointer_files=which(pointer_name,'-all');
if ischar(pointer_files)
    pointer_files={pointer_files};
end

pointer_file='';
megara_prefix=[megara_directory filesep];
for pointer_index=1:numel(pointer_files)
    candidate=pointer_files{pointer_index};
    if ~strncmp(candidate,megara_prefix,numel(megara_prefix))
        pointer_file=candidate;
        return
    end
end

end
