function remove_files(badfiles)

%Normalise input file list
badfiles=normalise_badfile_list(badfiles);

%Create removed files directory
removed_files_directory=fullfile(pwd,'removed_files');
if ~isfolder(removed_files_directory)
    mkdir(removed_files_directory);
end

%Set move path
movepath=removed_files_directory;

%move fit files to  folder
for file_index=1:numel(badfiles)
    file=badfiles{file_index};
    if isempty(file), continue; end
    if ~isfile(file), warning([file ' could not be moved because it does not exist.']); continue; end
    [move_status,move_message]=movefile(file,movepath,'f');
    if ~move_status
        disp([file ' could not be moved.'])
        warning(move_message)
        continue
    end
end

%Delete temporary g*.mat files
files=dir('g*.mat');
for file_index=1:numel(files)
    file_to_delete=fullfile(pwd,files(file_index).name);
    if isfile(file_to_delete), delete(file_to_delete); end
end

%remove mat files
files2=dir('bad*J*.mat');
for file_index=1:numel(files2)
    file_to_delete=fullfile(pwd,files2(file_index).name);
    if isfile(file_to_delete), delete(file_to_delete); end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function badfiles=normalise_badfile_list(badfiles)

%Normalise bad-file input
if isempty(badfiles)
    badfiles={};
elseif ischar(badfiles) || isstring(badfiles)
    badfiles=cellstr(string(badfiles));
elseif iscell(badfiles)
    badfiles=badfiles(:);
else
    error('badfiles must be a character vector, string, cell array, or empty.')
end

%Convert entries to character vectors
for file_index=1:numel(badfiles)
    if isstring(badfiles{file_index}), badfiles{file_index}=char(badfiles{file_index}); end
    if ~ischar(badfiles{file_index}), error('Each entry in badfiles must be a character vector or string scalar.'); end
end

end