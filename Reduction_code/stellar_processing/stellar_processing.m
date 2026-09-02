function stellar_processing(info,allorders,flatdata,jd_ths,fwave,m,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,objectname,past_red)

%Report processing start
disp(['processing stellar images @ ' datestr(now,'HH:MM:SS.FFF')])

%Record current raw-folder location
reduction_directory=pwd;

%Existing reductions in the current raw folder
existingMatFiles=dir(fullfile(reduction_directory,'J*.mat'));
existingMatNames={existingMatFiles.name};

%Keep only true reduced files: J1234567.mat, not J1234567_prc.mat
existingMatNames=existingMatNames(~cellfun(@isempty,regexp(existingMatNames,'^J\d+\.mat$','once')));
existingBaseNames=erase(existingMatNames,'.mat');

%Existing reductions already moved to the reduced target folder
if isempty(past_red)
    pastReducedNames={};
else
    pastReducedNames={past_red.name};
end

%Keep only true reduced files, not processing files
pastReducedNames=pastReducedNames(~cellfun(@isempty,regexp(pastReducedNames,'^J\d+\.mat$','once')));
pastReducedBaseNames=erase(pastReducedNames,'.mat');

%Only count target-folder reductions that correspond to files in this raw month/run
currentRawBaseNames=info.list(:);
pastReducedBaseNames=intersect(pastReducedBaseNames(:),currentRawBaseNames,'stable');

%Combined list of files that should not be reduced again for this raw month/run
alreadyReducedBaseNames=unique([existingBaseNames(:); pastReducedBaseNames(:)]);

%Report skipped observations as a current-folder summary only
fprintf('Found %d already reduced/rejected observations to skip for this raw folder\n',numel(alreadyReducedBaseNames))

%Build list of stellar files that are not already reduced
stellar_idx=[];
for file_info_index=1:numel(info.list)

    [exposure_type,jdexp]=make_temp_variable(info,file_info_index);

    %Skip non-stellar frames
    if ~strcmpi(exposure_type,'Stellar')
        continue
    end

    %Skip if already reduced in raw folder OR reduced target folder
    if any(strcmpi(info.list{file_info_index},alreadyReducedBaseNames))
        continue
    end

    %Skip if Julian date is missing
    if isempty(jdexp)
        continue
    end

    %Skip if this stellar frame is not the requested object
    if ~should_process_stellar_file(info,file_info_index,objectname)
        continue
    end

    %Keep this file for reduction
    stellar_idx(end+1)=file_info_index; %#ok<AGROW>

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMP TEST LIMIT - UNCOMMENT ONLY WHILE DEBUGGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%max_test_stellar_frames=3;
%stellar_idx=stellar_idx(1:min(max_test_stellar_frames,numel(stellar_idx)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMP SINGLE-FILE TEST - UNCOMMENT ONLY WHILE DEBUGGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test_stellar_file='J7850018';
%stellar_idx=stellar_idx(strcmp(info.list(stellar_idx),test_stellar_file));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Count stellar frames
nstellar=numel(stellar_idx);
fprintf('Found %d stellar frames to reduce\n',nstellar)

%Stop if no files need processing
if nstellar==0
    return
end

%Validate ThAr inputs
if isempty(jd_ths)
    error('No ThAr Julian dates available for stellar processing.')
end
if numel(jd_ths)<1
    error('At least one ThAr calibration is required for stellar processing.')
end

%Process stellar frames in a fresh pool for this raw folder
process_stellar_indices_in_parallel(stellar_idx,info,allorders,flatdata,jd_ths,fwave,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,reduction_directory)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function process_stellar_indices_in_parallel(stellar_idx,info,allorders,flatdata,jd_ths,fwave,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,reduction_directory)

%Process stellar frames in parallel using a fresh pool for this raw folder
number_stellar_files=numel(stellar_idx);

stellar_file_paths=cell(1,number_stellar_files);
stellar_display_names=cell(1,number_stellar_files);
jdexp_list=nan(1,number_stellar_files);
jd_th_before_list=nan(1,number_stellar_files);
jd_th_after_list=nan(1,number_stellar_files);
fwave1_list=cell(1,number_stellar_files);
fwave2_list=cell(1,number_stellar_files);

%Build absolute paths and ThAr calibration inputs before entering parfor
for stellar_index=1:number_stellar_files

    file_info_index=stellar_idx(stellar_index);
    [~,jdexp]=make_temp_variable(info,file_info_index);

    stellar_display_names{stellar_index}=[info.list{file_info_index} '.fit'];
    stellar_file_paths{stellar_index}=fullfile(reduction_directory,stellar_display_names{stellar_index});
    jdexp_list(stellar_index)=jdexp;

    [trueindbefore,trueindafter]=find_bracketing_thorium_indices(jd_ths,jdexp);

    if isempty(trueindbefore) && isempty(trueindafter)
        error(['No ThAr calibration images found for ' stellar_display_names{stellar_index}])
    end

    if isempty(trueindbefore)
        trueindbefore=trueindafter;
    end

    if isempty(trueindafter)
        trueindafter=trueindbefore;
    end

    [fwave1_list{stellar_index},fwave2_list{stellar_index}]=make_fwave_variable(trueindbefore,trueindafter,fwave);

    jd_th_before_list(stellar_index)=jd_ths(trueindbefore);
    jd_th_after_list(stellar_index)=jd_ths(trueindafter);

end

%Avoid creating a process pool and broadcasting the large flat structure
%when there is only one stellar frame in this observing folder.
if number_stellar_files==1
    fprintf('Stellar reductions: 0/1 complete\n')
    try
        run_stellar_frame_processing_in_folder(stellar_file_paths{1}, ...
            stellar_display_names{1},reduction_directory,info,allorders, ...
            flatdata,fwave1_list{1},fwave2_list{1},jd_th_before_list(1), ...
            jd_th_after_list(1),blue_data_chop_value,flatdata_medfilt, ...
            backdata_medfilt,cosmic_plotting,jdexp_list(1));
    catch ME
        error('Stellar processing failed for %s:\n%s', ...
            stellar_display_names{1}, ...
            getReport(ME,'extended','hyperlinks','off'))
    end
    fprintf('Stellar reductions: 1/1 complete\n')
    return
end

%Restart the process pool for this raw folder
parallel_cluster=parcluster('Processes');
number_parallel_workers=max(1,min(number_stellar_files, ...
    parallel_cluster.NumWorkers-1));
pool_object=restart_parallel_pool_for_reduction_folder(reduction_directory,number_parallel_workers);
pool_cleanup=onCleanup(@() delete(pool_object));

%Broadcast large read-only variables once per worker
info_constant=parallel.pool.Constant(info);
allorders_constant=parallel.pool.Constant(allorders);
flatdata_constant=parallel.pool.Constant(flatdata);

fprintf('Stellar reductions: 0/%d complete\n',number_stellar_files)
fprintf(['If no frame completes for 15 minutes, stop and rerun the reduction. ' ...
    'Completed frames will be skipped.\n'])

progress_count=0;
q=parallel.pool.DataQueue;
afterEach(q,@stellar_progress_update);

parfor stellar_index=1:number_stellar_files

    stellar_file_path=stellar_file_paths{stellar_index};
    stellar_display_name=stellar_display_names{stellar_index};
    jdexp=jdexp_list(stellar_index);
    jd_th_before=jd_th_before_list(stellar_index);
    jd_th_after=jd_th_after_list(stellar_index);
    fwave1=fwave1_list{stellar_index};
    fwave2=fwave2_list{stellar_index};

    try
        run_stellar_frame_processing_in_folder(stellar_file_path,stellar_display_name,reduction_directory,info_constant.Value,allorders_constant.Value,flatdata_constant.Value,fwave1,fwave2,jd_th_before,jd_th_after,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,jdexp);
    catch ME
        error('Stellar processing failed for %s:\n%s',stellar_display_name,getReport(ME,'extended','hyperlinks','off'))
    end

    send(q,stellar_display_name)

end

    function stellar_progress_update(stellar_display_name)

        %Update progress counter
        progress_count=progress_count+1;

        %Print compact progress only
        fprintf('Stellar reductions: %d/%d complete\n',progress_count,number_stellar_files)

    end

end

function run_stellar_frame_processing_in_folder(stellar_file_path,stellar_display_name,reduction_directory,info,allorders,flatdata,fwave1,fwave2,jd_th_before,jd_th_after,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,jdexp)

%Run stellar_frame_processing from the correct raw folder on the worker
if ~isfile(stellar_file_path)
    error(['Stellar FITS file does not exist: ' stellar_file_path])
end

%Keep legacy stellar_frame_processing filename behaviour while forcing worker folder state
original_worker_directory=pwd;
cleanup_object=onCleanup(@() cd(original_worker_directory));
cd(reduction_directory)

%Use the display name for the legacy processing function because it may build output filenames from the input string
stellar_frame_processing(stellar_display_name,info,allorders,flatdata,fwave1,fwave2,jd_th_before,jd_th_after,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,jdexp);

end

function pool_object=restart_parallel_pool_for_reduction_folder(reduction_directory,number_parallel_workers)

%Restart the process pool so workers do not retain state from a previous raw folder
existing_pool=gcp('nocreate');

if ~isempty(existing_pool)
    delete(existing_pool)
end

%Use the requested worker count
pool_object=parpool('Processes',number_parallel_workers);

%Process workers do not always inherit subfolders added to the client path
%after MATLAB starts. Add the bundled barycentric implementation explicitly
%before any stellar frame is processed.
stellar_processing_directory=fileparts(mfilename('fullpath'));
reduction_code_directory=fileparts(stellar_processing_directory);
matlab_barycentric_directory=fullfile(reduction_code_directory, ...
    'barycentric_correction','matlab_barycentric');
if ~isfolder(matlab_barycentric_directory)
    error('Bundled MATLAB barycentric-correction folder not found: %s', ...
        matlab_barycentric_directory)
end
addpath(matlab_barycentric_directory,'-begin')
path_futures=parfevalOnAll(@addpath,0,matlab_barycentric_directory,'-begin');
wait(path_futures)

%Set worker folders explicitly
folder_futures=parfevalOnAll(@cd,0,reduction_directory);
wait(folder_futures)

%Refresh worker function/path state
rehash_futures=parfevalOnAll(@rehash,0);
wait(rehash_futures)

end

function should_process=should_process_stellar_file(info,file_info_index,objectname)

%Accept all stellar files if no object filter is set
if isempty(objectname)
    should_process=true;
    return
end

%Use existing file information rather than rereading the FITS header
try
    [~,file_object_name]=identify_the_file(info,file_info_index);
catch
    warning(['Could not identify object name for ' info.list{file_info_index} ' from info structure. Skipping this file.'])
    should_process=false;
    return
end

%Compare requested and file object names
should_process=strcmpi(objectname,file_object_name);

end

function [trueindbefore,trueindafter]=find_bracketing_thorium_indices(jd_ths,jdexp)

%Find nearest ThAr before and after stellar exposure
jd_offsets=jd_ths(:)-jdexp;
trueindbefore=find(jd_offsets<=0,1,'last');
trueindafter=find(jd_offsets>=0,1,'first');

%Use nearest available side if exposure is outside ThAr range
if isempty(trueindbefore) && ~isempty(trueindafter)
    trueindbefore=trueindafter;
end

if isempty(trueindafter) && ~isempty(trueindbefore)
    trueindafter=trueindbefore;
end

end
