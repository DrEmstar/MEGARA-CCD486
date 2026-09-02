%BATCH_MOVE_FILES Move raw FITS files selected from file_info.mat metadata.
% This maintenance example moves Thorium frames longer than six seconds to
% removed_files. Run it from the raw-data directory, edit the selection
% test as needed, and review it carefully: movefile changes file locations.

clearvars
clc

load file_info.mat
if ~isfolder('removed_files')
    mkdir('removed_files')
end
for n=1:numel(info.list)
    file=info.list{n};
    exp_type=info.(file).HERCEXPT;
    exp_time=info.(file).EXPTIME;
    if strcmp(exp_type,'Thorium') && exp_time>6
        movefile([file '.fit'],'removed_files')
    end
end

