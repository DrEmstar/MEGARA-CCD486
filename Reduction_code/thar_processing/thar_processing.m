function [jd_ths,m,fwave,numbad_thar,save_flag,sum_info_rebuild]=thar_processing(info,allorders,Bias,flatdata,blue_data_chop_value,numbad_thar)

%ThAr processing
sum_info_rebuild=0;
save_flag=1;
jd_ths=[];
m=[];
fwave=struct();

%Minimum acceptable number of retained ThAr lines
minimum_kept_lines=800;

%Use absolute cache path so worker/current-folder state cannot affect this function
reduction_directory=pwd;
Thar_info_name=['ThAr_info_blue_' num2str(blue_data_chop_value) '.mat'];
Thar_info_path=fullfile(reduction_directory,Thar_info_name);

%Load cached ThAr information if available
if isfile(Thar_info_path)

    loaded_thar_information=load(Thar_info_path);

    required_variables={'jd_ths','m','fwave'};
    for required_variable_index=1:numel(required_variables)
        if ~isfield(loaded_thar_information,required_variables{required_variable_index})
            error(['Cached ThAr file ' Thar_info_name ' is missing variable ' required_variables{required_variable_index} '.'])
        end
    end

    %Check whether this is a MEGARA 2.0-compatible cache
    cache_has_version=isfield(loaded_thar_information,'megara_cache_version') || (isfield(loaded_thar_information,'fwave') && isstruct(loaded_thar_information.fwave) && isfield(loaded_thar_information.fwave,'megara_cache_version'));

    if cache_has_version

        jd_ths=loaded_thar_information.jd_ths;
        jd_ths=jd_ths(:).';
        m=loaded_thar_information.m;
        fwave=loaded_thar_information.fwave;
        save_flag=0;

        %Historical line counts cannot be recovered without rebuilding.
        if ~isfield(fwave,'quality_control')
            fwave.quality_control=struct();
            fwave.quality_control.cache_status='not_available_for_existing_cache';
            fwave.quality_control.messages={'Rebuild the ThAr cache to populate retained-line QC for existing frames.'};
            save_flag=1;
        end

        %Ensure version marker persists if the caller saves only fwave
        if ~isfield(fwave,'megara_cache_version')
            if isfield(loaded_thar_information,'megara_cache_version')
                fwave.megara_cache_version=loaded_thar_information.megara_cache_version;
            else
                fwave.megara_cache_version='MEGARA_2.0';
            end
        end

        %Process any ThAr frames that are present in the directory but missing from the cache
        thorium_info_indices=get_thorium_indices_from_info(info);
        existing_thar_count=count_cached_thars(fwave);
        directory_thar_count=numel(thorium_info_indices);

        if directory_thar_count>existing_thar_count

            missing_thorium_info_indices=thorium_info_indices(existing_thar_count+1:directory_thar_count);
            [new_jd_ths,new_m,new_fwave,numbad_thar_list,info_rebuild_list,new_quality_control]=process_thorium_indices_in_parallel(missing_thorium_info_indices,info,allorders,Bias,flatdata,blue_data_chop_value,minimum_kept_lines,reduction_directory);

            next_thar_number=existing_thar_count;

            for missing_thar_index=1:numel(missing_thorium_info_indices)

                if ~isfinite(new_jd_ths(missing_thar_index)) || isempty(new_fwave{missing_thar_index})
                    continue
                end

                next_thar_number=next_thar_number+1;
                jd_ths(next_thar_number)=new_jd_ths(missing_thar_index);
                fwave.(['th' num2str(next_thar_number)])=new_fwave{missing_thar_index};
                fwave.quality_control.(['th' num2str(next_thar_number)])=new_quality_control{missing_thar_index};

                if isempty(m) && ~isempty(new_m{missing_thar_index})
                    m=new_m{missing_thar_index};
                end

            end

            numbad_thar=numbad_thar+sum(numbad_thar_list);
            sum_info_rebuild=sum(info_rebuild_list);
            fwave.megara_cache_version='MEGARA_2.0';
            save_flag=1;

        end

        return

    else

        %Old cached ThAr files made before MEGARA 2.0.
        %Copy the old file aside, then rebuild a fresh cache.
        backup_Thar_info_path=create_legacy_thar_cache_path(reduction_directory,blue_data_chop_value);
        copyfile(Thar_info_path,backup_Thar_info_path,'f');

    end

end

%Process all ThAr files because no compatible cache was available
thorium_info_indices=get_thorium_indices_from_info(info);

if isempty(thorium_info_indices)
    warning('No ThAr frames found in this directory.')
    jd_ths=[];
    m=[];
    fwave=struct();
    fwave.megara_cache_version='MEGARA_2.0';
    save_flag=1;
    return
end

[jd_ths_t,m_t,fwave_t,numbad_thar_list,info_rebuild_list,quality_control_t]=process_thorium_indices_in_parallel(thorium_info_indices,info,allorders,Bias,flatdata,blue_data_chop_value,minimum_kept_lines,reduction_directory);

%Assemble ThAr outputs in directory order
valid_result_indices=find(isfinite(jd_ths_t) & ~cellfun(@isempty,fwave_t));

numbad_thar=numbad_thar+sum(numbad_thar_list);
sum_info_rebuild=sum(info_rebuild_list);

if isempty(valid_result_indices)
    warning('No ThAr frames were processed successfully.')
    jd_ths=[];
    m=[];
    fwave=struct();
    fwave.megara_cache_version='MEGARA_2.0';
    return
end

jd_ths=jd_ths_t(valid_result_indices);

for thar_number=1:numel(valid_result_indices)

    source_index=valid_result_indices(thar_number);

    fwave.(['th' num2str(thar_number)])=fwave_t{source_index};
    fwave.quality_control.(['th' num2str(thar_number)])=quality_control_t{source_index};

    if isempty(m) && ~isempty(m_t{source_index})
        m=m_t{source_index};
    end

end

%Mark new cache as MEGARA 2.0-compatible.
%This marker is inside fwave so it persists even if the caller saves only jd_ths, m, and fwave.
fwave.megara_cache_version='MEGARA_2.0';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [jd_ths_t,m_t,fwave_t,numbad_thar_list,info_rebuild_list,quality_control_t]=process_thorium_indices_in_parallel(thorium_info_indices,info,allorders,Bias,flatdata,blue_data_chop_value,minimum_kept_lines,reduction_directory)

%Process ThAr frames in parallel using a fresh pool for this raw folder
number_thorium_files=numel(thorium_info_indices);

jd_ths_t=nan(1,number_thorium_files);
numbad_thar_list=zeros(1,number_thorium_files);
info_rebuild_list=zeros(1,number_thorium_files);
fwave_t=cell(1,number_thorium_files);
m_t=cell(1,number_thorium_files);
quality_control_t=cell(1,number_thorium_files);

%Build absolute file paths before entering parfor
thar_file_paths=cell(1,number_thorium_files);
thar_display_names=cell(1,number_thorium_files);

for thorium_index=1:number_thorium_files
    file_info_index=thorium_info_indices(thorium_index);
    thar_display_names{thorium_index}=[info.list{file_info_index} '.fit'];
    thar_file_paths{thorium_index}=fullfile(reduction_directory,thar_display_names{thorium_index});
end

%Restart the process pool for this raw folder.
%This avoids stale worker folder/path state after the pipeline changes month.
parallel_cluster=parcluster('Processes');
number_parallel_workers=max(1,min(number_thorium_files, ...
    parallel_cluster.NumWorkers-1));
pool_object=restart_parallel_pool_for_reduction_folder(reduction_directory,number_parallel_workers);
pool_cleanup=onCleanup(@() delete(pool_object));

%Broadcast large read-only variables once per worker
info_constant=parallel.pool.Constant(info);
allorders_constant=parallel.pool.Constant(allorders);
flatdata_constant=parallel.pool.Constant(flatdata);
Bias_constant=parallel.pool.Constant(Bias);

fprintf('ThAr reductions: 0/%d complete\n',number_thorium_files)

progress_count=0;
q=parallel.pool.DataQueue;
afterEach(q,@thar_progress_update);

parfor thorium_index=1:number_thorium_files

    %Initialise temporary variables for parfor classification
    jd_ths_1=NaN;
    m_temp=[];
    fwave_temp=[];
    numbad_thar_flag=0;
    info_rebuild=0;
    kept_line_count=NaN;

    thar_file_path=thar_file_paths{thorium_index};
    thar_display_name=thar_display_names{thorium_index};

    try
        [jd_ths_1,m_temp,fwave_temp,numbad_thar_flag,info_rebuild,kept_line_count]=process_thar_frames(thar_file_path,info_constant.Value,allorders_constant.Value,Bias_constant.Value,flatdata_constant.Value,blue_data_chop_value);
    catch ME
        error('ThAr processing failed for %s:\n%s',thar_display_name,getReport(ME,'extended','hyperlinks','off'))
    end

    jd_ths_t(thorium_index)=jd_ths_1;
    numbad_thar_list(thorium_index)=numbad_thar_flag;
    fwave_t{thorium_index}=fwave_temp;
    m_t{thorium_index}=m_temp;
    info_rebuild_list(thorium_index)=info_rebuild;

    %Save the line-count quality control already produced by calibration.
    thar_quality_control=struct();
    thar_quality_control.version='MEGARA_2.1_QC_1';
    thar_quality_control.stage='thorium_wavelength_calibration';
    thar_quality_control.file_name=thar_display_name;
    thar_quality_control.status='pass';
    thar_quality_control.messages={};
    thar_quality_control.kept_line_count=kept_line_count;
    thar_quality_control.minimum_recommended_line_count=minimum_kept_lines;
    if isnan(kept_line_count)
        thar_quality_control.status='warning';
        thar_quality_control.messages{end+1}='The retained ThAr line count could not be measured.';
    elseif kept_line_count<minimum_kept_lines
        thar_quality_control.status='warning';
        thar_quality_control.messages{end+1}=sprintf('Only %.0f ThAr lines were retained; the current recommendation is at least %.0f.',kept_line_count,minimum_kept_lines);
    end
    quality_control_t{thorium_index}=thar_quality_control;

    %Send simple cell message only; avoid struct assignment inside parfor
    send(q,{thar_display_name,kept_line_count})

end

    function thar_progress_update(progress_message)

        %Update progress counter
        progress_count=progress_count+1;

        %Extract message contents
        file_name=progress_message{1};
        kept_line_count=progress_message{2};

        %Warn only for low retained-line count
        if ~isnan(kept_line_count) && kept_line_count<minimum_kept_lines
            fprintf(2,'WARNING: %s kept only %.0f ThAr lines.\n',file_name,kept_line_count)
        end

        %Print compact progress only
        fprintf('ThAr reductions: %d/%d complete\n',progress_count,number_thorium_files)

    end

end

function pool_object=restart_parallel_pool_for_reduction_folder(reduction_directory,number_parallel_workers)

%Restart the process pool so workers do not retain state from a previous raw folder
existing_pool=gcp('nocreate');

if ~isempty(existing_pool)
    delete(existing_pool)
end

%Use the requested worker count
pool_object=parpool('Processes',number_parallel_workers);

%Set worker folders explicitly
folder_futures=parfevalOnAll(@cd,0,reduction_directory);
wait(folder_futures)

%Refresh worker function/path state
rehash_futures=parfevalOnAll(@rehash,0);
wait(rehash_futures)

end

function backup_Thar_info_path=create_legacy_thar_cache_path(reduction_directory,blue_data_chop_value)

%Create legacy ThAr-cache backup path
backup_Thar_info_name=['ThAr_info_blue_' num2str(blue_data_chop_value) '_Megara_1.6.mat'];
backup_Thar_info_path=fullfile(reduction_directory,backup_Thar_info_name);

%Avoid overwriting an existing backup
if isfile(backup_Thar_info_path)
    timestamp=datestr(now,'yyyymmdd_HHMMSS');
    backup_Thar_info_name=['ThAr_info_blue_' num2str(blue_data_chop_value) '_Megara_1.6_' timestamp '.mat'];
    backup_Thar_info_path=fullfile(reduction_directory,backup_Thar_info_name);
end

end

function cached_thar_count=count_cached_thars(fwave)

%Count cached fwave.thN fields
cached_thar_count=0;

while true

    candidate_field_name=['th' num2str(cached_thar_count+1)];

    if isfield(fwave,candidate_field_name)
        cached_thar_count=cached_thar_count+1;
    else
        break
    end

end

end

function thorium_info_indices=get_thorium_indices_from_info(info)

%Find ThAr frames in info list
thorium_info_indices=[];

for file_info_index=1:numel(info.list)

    [exposure_type,~,~]=make_temp_variable_thars(info,file_info_index);

    if strcmpi(exposure_type,'Thorium')
        thorium_info_indices(end+1)=file_info_index; %#ok<AGROW>
    end

end

end
