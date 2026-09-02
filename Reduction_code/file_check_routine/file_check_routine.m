function [change_flag,counter_good]=file_check_routine( ...
    starname,skip_manual_file_check)

if nargin<2 || isempty(skip_manual_file_check)
    skip_manual_file_check=false;
end

%%%%%%   PARAMETERS   %%%%%%%
%leave these, used later
badfiles=[];
readout=[];

%%%%%%   MAIN SECTION %%%%%%%

if skip_manual_file_check
    disp(['Manual Good/Bad file review is disabled; ambiguous readable ' ...
        'files will be retained.'])
end

%general bad file routine
has_run=badfile_check_general;

%Load current file information
if isfile('file_info.mat')
    loaded_file_info=load('file_info.mat','info');
    if isfield(loaded_file_info,'info'), info=loaded_file_info.info; else, error('file_info.mat exists but does not contain variable info.'); end
else
    info=display_directory_file_information_nooutput();
    save('file_info','info')
end

%Track checks separately from changes to the raw directory. Rebuilding a
%missing classification does not require all FITS headers to be reread.
check_performed=false;
change_flag=false;

%CYCLE THROUGH FILES FIRST CHECKS FOR CALIBRATIONS AND SKIPS IF NOT
%Get line profile and move it to the right folder
test_badcal=dir('badfiles_calibration.mat');
if numel(test_badcal)>0
    disp('Calibration file check already completed in this directory')
    badfiles_calibration=[];
else
    check_performed=true;
    waitbar_handle=waitbar(0,'Please wait...');
    cleanup_waitbar=onCleanup(@() close_waitbar_if_open(waitbar_handle));
    counter_good=0;
    counter_bad=0;
    counter_more_good=0;
    previous_good=0;
    othergoodtest.v1=[];
    othergoodtype.v1=[];
    otherbadtest.v1=[];
    disp('Testing calibration files')
    badfiles_calibration=[];
    for file_info_index=1:numel(info.list)
        waitbar(file_info_index/max(numel(info.list),1),waitbar_handle,info.list{file_info_index})
        %I.D. the file
        [type,star]=identify_the_file(info,file_info_index);
        if has_existing_line_analysis(type,info.list{file_info_index})
            continue
        end
        if strcmp(type,'White L') || strcmp(type,'Thorium')
            try
                test_line=read_line_profile_with_retry( ...
                    star,starname,type,file_info_index,info);
            catch profile_error
                [badfiles_calibration,readout,counter_bad]= ...
                    record_unreadable_fits(type,info,file_info_index, ...
                    profile_error,badfiles_calibration,readout,counter_bad);
                continue
            end
            clear good
            [badfiles_calibration,readout,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good]=analyse_line_profile(type,test_line,badfiles_calibration,readout,info,file_info_index,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good,skip_manual_file_check);
        end
    end
    char(readout);
    close_waitbar_if_open(waitbar_handle)
    save('badfiles_calibration','badfiles_calibration')
end

%CYCLE THROUGH FILES FIRST CHECKS FOR STELLARS AND SKIPS IF NOT
%Get line profile and move it to the right folder
counter_good=0;
counter_more_good=0;
previous_good=0;
counter_bad=0;
all_badfiles_stellar=[];
%include alias names
alt_starnames=get_alternative_star_names(starname);

%Check stellar files for each alias
for alternative_name_index=1:numel(alt_starnames)
    current_starname=alt_starnames{alternative_name_index};
    badstellar_name=['badfiles_' current_starname '.mat'];
    test_badstellar=dir(badstellar_name);

    %A legacy completion marker is not sufficient if its per-frame
    %classifications have since been removed or were never saved.
    if ~isempty(test_badstellar) && ...
            ~stellar_check_cache_is_complete(info,current_starname)
        disp(['Cached stellar file check for ' current_starname ...
            ' is incomplete; rebuilding the missing classifications'])
        test_badstellar=[];
    end

    if numel(test_badstellar)>0
        disp(['Stellar file check ' current_starname ' already completed in this directory'])
        for file_info_index=1:numel(info.list)
            %I.D. the file
            [type,star]=identify_the_file(info,file_info_index);
            raw_file_name=[info.list{file_info_index} '.fit'];
            if strcmp(type,'Stellar') && strcmpi(star,current_starname) && ...
                    existing_line_analysis_is_good(type, ...
                    info.list{file_info_index}) && isfile(raw_file_name)
                counter_good=counter_good+1;
            end
        end
    else
        check_performed=true;
        waitbar_handle=waitbar(0,'Please wait...');
        cleanup_waitbar=onCleanup(@() close_waitbar_if_open(waitbar_handle));
        badfiles_stellar=[];
        othergoodtest.v1=[];
        othergoodtype.v1=[];
        otherbadtest.v1=[];
        for file_info_index=1:numel(info.list)
            waitbar(file_info_index/max(numel(info.list),1),waitbar_handle,info.list{file_info_index})
            %I.D. the file
            [type,star]=identify_the_file(info,file_info_index);
            if has_existing_line_analysis(type,info.list{file_info_index})
                if strcmp(type,'Stellar') && strcmpi(star,current_starname) && ...
                        existing_line_analysis_is_good(type, ...
                        info.list{file_info_index}) && ...
                        isfile([info.list{file_info_index} '.fit'])
                    counter_good=counter_good+1;
                end
                continue
            end
            if strcmp(type,'Stellar') && strcmpi(star,current_starname) && ...
                    isfile([info.list{file_info_index} '.fit'])
                try
                    test_line=read_line_profile_with_retry( ...
                        star,current_starname,type,file_info_index,info);
                catch profile_error
                    [badfiles_stellar,readout,counter_bad]= ...
                        record_unreadable_fits(type,info,file_info_index, ...
                        profile_error,badfiles_stellar,readout,counter_bad);
                    continue
                end
                clear good
                [badfiles_stellar,readout,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good]=analyse_line_profile(type,test_line,badfiles_stellar,readout,info,file_info_index,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good,skip_manual_file_check);
            end
        end
        char(readout);
        close_waitbar_if_open(waitbar_handle)
        save(badstellar_name,'badfiles_stellar')
        all_badfiles_stellar=cat(1,all_badfiles_stellar,badfiles_stellar);
    end
end

%Add extra accepted files
counter_good=counter_good+counter_more_good;

%%%%%%%%%%%%%%%%%%%  MAIN REDUCTION COMPLETE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE FOLLOWING WILL PLOT EACH OF THE SUB FOLDERS IN LINE ANALYSIS %%%%%%
%%%%%%%%%%%%%%%%%%%%%% SAVE THEM. THIS IS FOR YOU TO INSPECT %%%%%%%%%%%%%%

%Find line-analysis files
god_names=dir('good_*.mat');
bad_names=dir('bad_*.mat');

%Plot accepted and rejected line profiles
if check_performed && ~skip_manual_file_check && ...
        (isempty(bad_names)==0 || isempty(god_names)==0)
    figure(1000)
    for file_index=1:numel(bad_names)
        loaded_bad_file=load(bad_names(file_index).name);
        if ~isfield(loaded_bad_file,'test_line'), continue; end
        classification_file_name=bad_names(file_index).name;
        if numel(classification_file_name)>=5 && classification_file_name(5)=='S'
            subplot(2,3,4)
            plot(loaded_bad_file.test_line)
            hold on
        elseif numel(classification_file_name)>=5 && classification_file_name(5)=='W'
            subplot(2,3,5)
            plot(loaded_bad_file.test_line)
            hold on
        elseif numel(classification_file_name)>=5 && classification_file_name(5)=='T'
            subplot(2,3,6)
            plot(loaded_bad_file.test_line)
            hold on
        end
    end

    for file_index=1:numel(god_names)
        loaded_good_file=load(god_names(file_index).name);
        if ~isfield(loaded_good_file,'test_line'), continue; end
        classification_file_name=god_names(file_index).name;
        if numel(classification_file_name)>=6 && classification_file_name(6)=='S'
            subplot(2,3,1)
            plot(loaded_good_file.test_line)
            hold on
        elseif numel(classification_file_name)>=6 && classification_file_name(6)=='W'
            subplot(2,3,2)
            plot(loaded_good_file.test_line)
            hold on
        elseif numel(classification_file_name)>=6 && classification_file_name(6)=='T'
            subplot(2,3,3)
            plot(loaded_good_file.test_line)
            hold on
        end
    end

    subplot(2,3,1)
    xlabel(strcat(num2str(numel(dir('good_S*.mat'))),' ',' -Stellar Accepted'))
    subplot(2,3,2)
    xlabel(strcat(num2str(numel(dir('good_W*.mat'))),' ','-White Lamps Accepted'))
    subplot(2,3,3)
    xlabel(strcat(num2str(numel(dir('good_T*.mat'))),' ','-Thoriums Accepted'))
    subplot(2,3,4)
    xlabel(strcat(num2str(numel(dir('bad_S*.mat'))),' ','-Rejected'))
    subplot(2,3,5)
    xlabel(strcat(num2str(numel(dir('bad_W*.mat'))),' ','-Rejected'))
    subplot(2,3,6)
    xlabel(strcat(num2str(numel(dir('bad_T*.mat'))),' ','-Rejected'))

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Press enter to move bad files to a subdirectory or CTRL-C to change')
    user_input=input('break out.');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    close_figure_if_open(1000)
elseif check_performed && skip_manual_file_check
    disp(['Manual Good/Bad review skipped. Automatically rejected files ' ...
        'will be moved to removed_files.'])
else
end

%Only moving FITS files changes the directory information cache.
change_flag=~isempty(badfiles_calibration) || ...
    ~isempty(all_badfiles_stellar);

%Move rejected calibration files
if exist('badfiles_calibration','var') && ~isempty(badfiles_calibration)
    remove_files(badfiles_calibration)
end

%Move rejected stellar files
if exist('all_badfiles_stellar','var') && ~isempty(all_badfiles_stellar)
    remove_files(all_badfiles_stellar)
end

%Report completion
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('The file check program is now complete. The bad files have been')
disp('moved to a separate folder called removed_files and should not be reduced.')
disp('Everything else remaining in the current directory is suitable for further')
disp('use. Reduction will commence shortly.')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function analysis_exists=has_existing_line_analysis(type,file_identifier)

%Check for previous line-analysis file
if isempty(type) || isempty(file_identifier)
    analysis_exists=false;
    return
end
analysis_exists=false;
legacy_file_name=['god' type file_identifier '.mat'];
if isfile(legacy_file_name), analysis_exists=true; return; end
good_matches=dir(['good_*' file_identifier '*.mat']);
bad_matches=dir(['bad_*' file_identifier '*.mat']);
if ~isempty(good_matches) || ~isempty(bad_matches), analysis_exists=true; end

end

function analysis_is_good=existing_line_analysis_is_good(type,file_identifier)

analysis_is_good=false;
legacy_file_name=['god' type file_identifier '.mat'];
if isfile(legacy_file_name)
    analysis_is_good=true;
    return
end

good_matches=dir(['good_*' file_identifier '*.mat']);
if ~isempty(good_matches)
    analysis_is_good=true;
end

end

function cache_is_complete=stellar_check_cache_is_complete(info,starname)

%Every matching raw stellar frame still present must have a classification
cache_is_complete=true;
for file_info_index=1:numel(info.list)
    [type,star]=identify_the_file(info,file_info_index);
    raw_file_name=[info.list{file_info_index} '.fit'];
    if strcmp(type,'Stellar') && strcmpi(star,starname) && ...
            isfile(raw_file_name) && ...
            ~has_existing_line_analysis(type,info.list{file_info_index})
        cache_is_complete=false;
        return
    end
end

end

function alt_starnames=get_alternative_star_names(starname)

%Initialise aliases
alt_starnames={starname};

%Load known aliases
if ~isfile('Mus_known_alias_names.mat')
    return
end
loaded_aliases=load('Mus_known_alias_names.mat');
if ~isfield(loaded_aliases,'wildc')
    return
end
wildc=loaded_aliases.wildc;

%Find aliases
alias_rows=find(any(strcmpi(wildc,starname),2));
for alias_row_index=1:numel(alias_rows)
    alias_row=alias_rows(alias_row_index);
    for alias_column_index=1:size(wildc,2)
        alias_name=wildc{alias_row,alias_column_index};
        if isempty(alias_name), continue; end
        if ~any(strcmpi(alt_starnames,alias_name))
            alt_starnames{end+1}=alias_name;
        end
    end
end

end

function close_waitbar_if_open(waitbar_handle)

%Close waitbar safely
if exist('waitbar_handle','var') && ishghandle(waitbar_handle)
    close(waitbar_handle)
end

end

function close_figure_if_open(figure_number)

%Close figure safely
figure_handle=findobj('Type','figure','Number',figure_number);
if ~isempty(figure_handle)
    close(figure_handle)
end

end

function test_line=read_line_profile_with_retry( ...
    star,starname,type,file_info_index,info)

%Retry once in case an external drive has a temporary read failure.
try
    test_line=get_the_line_profile( ...
        star,starname,type,file_info_index,info);
catch
    pause(0.5)
    test_line=get_the_line_profile( ...
        star,starname,type,file_info_index,info);
end

end

function [badfiles,readout,counter_bad]=record_unreadable_fits( ...
    type,info,file_info_index,profile_error,badfiles,readout,counter_bad)

%Record an unreadable image as bad without stopping the complete target.
file_root=info.list{file_info_index};
file_name=[file_root '.fit'];
file_path=fullfile(pwd,file_name);
warning('MEGARA:UnreadableFITS', ...
    ['Could not read image data from ' file_name ...
    ' after two attempts. The file will be treated as bad.\n' ...
    profile_error.message])
badfiles=cat(1,badfiles,{file_path});
readout=cat(1,readout,{[file_root ' a ' type ' is unreadable']});
counter_bad=counter_bad+1;
test_line=NaN;
read_error_message=profile_error.message;
safe_type=regexprep(type,'\s+','_');
save(['bad_' safe_type file_root '.mat'], ...
    'test_line','read_error_message')

end
