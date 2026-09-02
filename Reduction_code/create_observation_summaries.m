function [passback_objectname]=create_observation_summaries(raw_data_location,reduced_data_location,reduction_folders,reduction_objectname)

%Initialise census settings
build_database=false;%If this is turned on the database is built from all stellar files in the directory
wildcard_matching=false;
database_build_counter=0;
database_store=struct();
passback_objectname=reduction_objectname;

%Store current directory
starting_directory=pwd;
cleanup_directory=onCleanup(@() cd(starting_directory));

%Load wildcard aliases
wildc=load_wildcard_aliases();

%Build target list
[database_obs,wildcard_matching,build_database,passback_objectname]=initialise_database_targets(reduction_objectname);

%Initialise counters
port_4_count=zeros(max(numel(database_obs),1000),1);
fibre_3_count=zeros(max(numel(database_obs),1000),1);

% Get a list of all observing folders
subFolders=get_selected_observing_folders(raw_data_location,reduction_folders);

% Collect all the file_information in all folders
port_flag=0;
fibre_flag=0;

%Scan observing folders
for folder_index=1:numel(subFolders)
    fprintf('Running census on Month/Run %.0f of %.0f\n',folder_index,numel(subFolders))
    current_raw_directory=fullfile(raw_data_location,subFolders(folder_index).name);
    cd(current_raw_directory)
    fibre_error_num=0;
    fibre_error_list=struct();

    %Load or create file information
    if isfile('file_info.mat')
        loaded_file_info=load('file_info.mat');
        if isfield(loaded_file_info,'info'), info=loaded_file_info.info; else, error(['file_info.mat in ' current_raw_directory ' does not contain variable info.']); end
    else
        info=display_directory_file_information_nooutput();
        save('file_info','info')
    end

    %skip checking again if reduction already happened here
    if isfolder('removed_files') || isfolder('Removed files')
    else
        %general bad file check
        badfile_check_general;
        file_list=dir('J*.fit');

        %Move unsupported files out of the active raw folder
        for file_index=1:numel(file_list)
            header=get_header(file_list(file_index).name);
            temp_head=get_required_headers_from_header(header);
            port_value=temp_head.PORT;
            fibre_value=temp_head.HERCFIB;
            if ~isempty(port_value) && port_value>1
                port_flag=port_flag+1;
                if port_flag==1 && ~isfolder('four_port_files'), mkdir('four_port_files'); end
                movefile(file_list(file_index).name,'four_port_files')
            end
            if ~isempty(fibre_value) && fibre_value==3
                fibre_flag=fibre_flag+1;
                if fibre_flag==1 && ~isfolder('fibre_3_files'), mkdir('fibre_3_files'); end
                movefile(file_list(file_index).name,'fibre_3_files')
            elseif ~isempty(fibre_value) && fibre_value==1
                continue
                %leave fibre error files in dir
            else
                fibre_error_num=fibre_error_num+1;
                fibre_error_list.(['name' num2str(fibre_error_num)])=file_list(file_index).name;
            end
        end
        info=display_directory_file_information_nooutput();
        save('file_info','info')
    end

    %Skip empty folders
    if ~isfield(info,'list') || isempty(info.list), continue; end

    %Track which targets appear in this folder
    target_seen_this_folder=false(max(numel(database_obs),1000),1);

    %Collect observations for selected targets
    for file_info_index=1:length(info.list)
        [type,star]=identify_the_file(info,file_info_index); %I.D. the file
        star=normalise_object_name(star);

        %Build database from all stellar targets if requested
        if build_database==1
            if strcmpi(type,'Stellar') && ~any(strcmpi(database_obs,star))
                database_build_counter=database_build_counter+1;
                database_obs{database_build_counter,1}=star;
                passback_objectname.(['target' num2str(database_build_counter)])=star;
                if database_build_counter>numel(port_4_count), port_4_count(database_build_counter,1)=0; fibre_3_count(database_build_counter,1)=0; target_seen_this_folder(database_build_counter,1)=false; end
            end
        end

        %Update each target database
        for target_index=1:length(database_obs)
            objectname=normalise_object_name(database_obs{target_index});
            if isempty(objectname), continue; end
            objectname2=get_hd_alias_name(objectname);
            if strcmpi(type,'Stellar') && (strcmpi(star,objectname) || strcmpi(star,objectname2))
                observation_entry=get_observation_entry(info,file_info_index);
                database_store=append_observation_to_database(database_store,objectname,observation_entry,subFolders(folder_index).name,target_seen_this_folder(target_index));
                target_seen_this_folder(target_index)=true;
            end
        end
    end
end

%wildcard matching for musician targets
if wildcard_matching==1
    database_store=merge_wildcard_alias_databases(database_store,wildc,database_obs);
end

% if port_flag>0
%     fprintf('A total of %.0f observations are 4-port readout and not reduced\n',port_flag)
% end
% if fibre_flag>0
%     fprintf('A total of %.0f observations are fibre 3 and not reduced\n',fibre_flag)
% end
% if fibre_error > 0;
%     fprintf(' A total of %.0f observations have no recorded fibre position and may not reduce\n',fibre_error)
% end

%Create observation summary directory
observation_summary_directory=fullfile(reduced_data_location,'Observation Summaries');
if ~isfolder(observation_summary_directory), mkdir(observation_summary_directory); end
addpath(observation_summary_directory)
cd(observation_summary_directory)

%Write observation summary files
master_file=cell(length(database_obs),4);
for target_index=1:length(database_obs)
    objectname=normalise_object_name(database_obs{target_index});
    master_file{target_index,1}=objectname;
    database_key=get_database_key(objectname);
    if isfield(database_store,database_key)
        star_database=database_store.(database_key);
        if isfield(star_database,'months'), star_database.months=unique(star_database.months,'stable'); end
        master_file{target_index,2}=numel(star_database.filename);
        master_file{target_index,3}=port_4_count(min(target_index,numel(port_4_count)));
        master_file{target_index,4}=fibre_3_count(min(target_index,numel(fibre_3_count)));
        save_star_database(objectname,star_database)
    end
end

%Save master summary
save('observation_summary_count','master_file')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wildc=load_wildcard_aliases()

%Load known aliases
wildc={};
if isfile('Mus_known_alias_names.mat')
    loaded_aliases=load('Mus_known_alias_names.mat');
    if isfield(loaded_aliases,'wildc'), wildc=loaded_aliases.wildc; end
end

end

function [database_obs,wildcard_matching,build_database,passback_objectname]=initialise_database_targets(reduction_objectname)

%Initialise target mode
database_obs={};
wildcard_matching=false;
build_database=false;
passback_objectname=reduction_objectname;

%Load predefined target database
if ischar(reduction_objectname) || isstring(reduction_objectname)
    reduction_objectname=char(reduction_objectname);
    if strcmp(reduction_objectname,'musician') || strcmp(reduction_objectname,'tess') || strcmp(reduction_objectname,'exocomet')
        if strcmp(reduction_objectname,'musician')
            loaded_database=load('mus_database_cool.mat','mus_database');
            database_obs=loaded_database.mus_database(:,1);
        end
        if strcmp(reduction_objectname,'tess')
            loaded_database=load('tess_database.mat','tess_database');
            database_obs=loaded_database.tess_database(:,1);
        end
        if strcmp(reduction_objectname,'exocomet')
            loaded_database=load('exocomet_database.mat','exocomet_database');
            database_obs=loaded_database.exocomet_database(:,1);
        end
        wildcard_matching=true;
    elseif isempty(reduction_objectname)
        database_obs={};
        build_database=true;
        passback_objectname=struct();
    end
elseif isa(reduction_objectname,'struct')
    target_field_names=fieldnames(reduction_objectname);
    database_obs=cell(numel(target_field_names),1);
    for target_index=1:numel(target_field_names)
        expected_target_field_name=['target' num2str(target_index)];
        if isfield(reduction_objectname,expected_target_field_name)
            database_obs{target_index}=reduction_objectname.(expected_target_field_name);
        else
            database_obs{target_index}=reduction_objectname.(target_field_names{target_index});
        end
    end
elseif isempty(reduction_objectname)
    database_obs={};
    build_database=true;
    passback_objectname=struct();
else
    error('reduction_objectname must be a character vector, string scalar, structure, or empty.')
end

%Normalise target names
for target_index=1:numel(database_obs)
    database_obs{target_index}=normalise_object_name(database_obs{target_index});
end

end

function subFolders=get_selected_observing_folders(raw_data_location,reduction_folders)

%Read observing folders
files=dir(raw_data_location);
directory_flags=[files.isdir];
subFolders=files(directory_flags);
subFolders=subFolders(~ismember({subFolders.name},{'.','..'}));

%set up reduction folder subsets
if numel(reduction_folders)>0
    selected_subFolders=struct('name',{},'folder',{},'date',{},'bytes',{},'isdir',{},'datenum',{});
    folder_names={subFolders.name};
    if ischar(reduction_folders) || isstring(reduction_folders), reduction_folders=cellstr(reduction_folders); end
    for folder_index=1:numel(reduction_folders)
        folder_match=strcmp(reduction_folders{folder_index},folder_names);
        hits=find(folder_match);
        if numel(hits)==0
            disp(['Folder ' reduction_folders{folder_index} ' does not exist in the raw data directory, please check name and directory location'])
        else
            selected_subFolders(end+1)=subFolders(hits(1));
        end
    end
    subFolders=selected_subFolders;
end

end

function objectname=normalise_object_name(objectname)

%Normalise target name
if isstring(objectname), objectname=char(objectname); end
if isempty(objectname), return; end
objectname=strrep(objectname,'-','_');

end

function objectname2=get_hd_alias_name(objectname)

%Create HD alias
if numel(objectname)>=4
    objectname2=['HD' objectname(4:end)];
else
    objectname2=objectname;
end

end

function observation_entry=get_observation_entry(info,file_info_index)

%Extract observation metadata
file_root=info.list{file_info_index};
file_information=info.(file_root);
observation_entry.filename=[file_root '.fit'];
observation_entry.date=get_optional_field(file_information,'DATE_OBS');
observation_entry.jd=get_optional_field(file_information,'MJD_OBS');
observation_entry.exptime=get_optional_field(file_information,'EXPTIME');

end

function value=get_optional_field(input_structure,field_name)

%Read optional field
if isfield(input_structure,field_name)
    value=input_structure.(field_name);
else
    value=[];
end

end

function database_store=append_observation_to_database(database_store,objectname,observation_entry,month_name,target_already_seen_this_folder)

%Append observation to database
database_key=get_database_key(objectname);
if ~isfield(database_store,database_key), database_store.(database_key)=initialise_star_database(); end
star_database=database_store.(database_key);
next_observation_index=numel(star_database.filename)+1;
star_database.filename{next_observation_index}=observation_entry.filename;
star_database.date{next_observation_index}=observation_entry.date;
star_database.jd{next_observation_index}=observation_entry.jd;
star_database.exptime(next_observation_index)=observation_entry.exptime;
if ~target_already_seen_this_folder
    star_database.months{numel(star_database.months)+1}=month_name;
end
database_store.(database_key)=star_database;

end

function star_database=initialise_star_database()

%Create empty star database
star_database.filename={};
star_database.date={};
star_database.jd={};
star_database.exptime=[];
star_database.months={};

end

function database_key=get_database_key(objectname)

%Create valid database key
database_key=matlab.lang.makeValidName(['star_database_' objectname]);

end

function database_store=merge_wildcard_alias_databases(database_store,wildc,database_obs)

%Merge wildcard alias databases
if isempty(wildc), return; end
for wildcard_index=1:size(wildc,1)
    wildn=normalise_object_name(wildc{wildcard_index,1});
    wildn_key=get_database_key(wildn);
    if ~isfield(database_store,wildn_key), database_store.(wildn_key)=initialise_star_database(); end
    for alias_index=2:size(wildc,2)
        if isempty(wildc{wildcard_index,alias_index}), continue; end
        wilda=normalise_object_name(wildc{wildcard_index,alias_index});
        wilda_key=get_database_key(wilda);
        if isfield(database_store,wilda_key)
            database_store.(wildn_key)=append_star_database(database_store.(wildn_key),database_store.(wilda_key));
            if ~strcmpi(wildn_key,wilda_key), database_store=rmfield(database_store,wilda_key); end
        end
    end
end

end

function target_database=append_star_database(target_database,alias_database)

%Append alias database
if isempty(alias_database.filename), return; end
start_index=numel(target_database.filename);
for observation_index=1:numel(alias_database.filename)
    target_database.filename{start_index+observation_index}=alias_database.filename{observation_index};
    target_database.date{start_index+observation_index}=alias_database.date{observation_index};
    target_database.jd{start_index+observation_index}=alias_database.jd{observation_index};
    target_database.exptime(start_index+observation_index)=alias_database.exptime(observation_index);
end
target_database.months=unique([target_database.months alias_database.months],'stable');

end

function save_star_database(objectname,star_database)

%Save database using legacy variable name
database_variable_name=matlab.lang.makeValidName(['star_database_' objectname]);
database_file_name=['star_database_' objectname '.mat'];
database_output=struct();
database_output.(database_variable_name)=star_database;
save(database_file_name,'-struct','database_output')

end