function has_run=badfile_check_general()
%see if this has already run- skip if not

%Initialise outputs
has_run=false;
badfiles_general={};
readout1={};
readout2={};

%Skip if check has already run
existing_check_file=dir('badfiles_general.mat');
if numel(existing_check_file)>0
    %disp('General file check already completed in this directory')
    has_run=true;
    return
end

%remove duplicate files sometimes generated
duplicate_file_list=dir('J*@*.fit');
if ~isempty(duplicate_file_list)
    for duplicate_file_index=1:length(duplicate_file_list)
        if isfile(duplicate_file_list(duplicate_file_index).name), delete(duplicate_file_list(duplicate_file_index).name); end
    end
end

%Set bad-file size limits
minimum_valid_file_size_bytes=33842880;
maximum_valid_file_size_bytes=34000000;

%Check FITS file sizes
file_size_list=dir('J*.fit');

%this checks file size. anything <33.8mb will be bad!!!!
for file_index=1:numel(file_size_list)
    current_file_name=file_size_list(file_index).name;
    current_file_path=fullfile(pwd,current_file_name);
    if file_size_list(file_index).bytes<minimum_valid_file_size_bytes
        badfiles_general=cat(1,badfiles_general,{current_file_path});
        readout1=cat(1,readout1,{[current_file_name ' is an incomplete file.']});
        continue
    end
    if file_size_list(file_index).bytes>maximum_valid_file_size_bytes
        badfiles_general=cat(1,badfiles_general,{current_file_path});
        readout2=cat(1,readout2,{[current_file_name ' is a large file.']});
        continue
    end
end

%Move bad files
if ~isempty(badfiles_general)
    remove_files(badfiles_general)
end

%Refresh file information
info=display_directory_file_information_nooutput();
save('file_info','info')

%Save check record
save('badfiles_general','badfiles_general','readout1','readout2','minimum_valid_file_size_bytes','maximum_valid_file_size_bytes')
has_run=true;

end