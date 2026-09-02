function [badfiles,readout,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good]=analyse_line_profile(type,test_line,badfiles,readout,files,s,counter_good,counter_bad,othergoodtest,othergoodtype,otherbadtest,counter_more_good,previous_good,skip_manual_file_check)

if nargin<14 || isempty(skip_manual_file_check)
    skip_manual_file_check=false;
end

%Measure profile strength
test=mean(test_line,'omitnan');

%Initialise decision state
good=[];
think=[];
manual_override=false;

%Initial boundaries for defining good files.
if test>0.014 && test<0.025
    think='Thorium';
    good=strcmp(type,think);
elseif test>0.035 && test<0.04
    think='White L';
    good=strcmp(type,think);
elseif test>0.065 && test<0.08
    think='Stellar';
    good=strcmp(type,think);
    if good, counter_good=counter_good+1; end
else
    %Here are tests based on previous good identifications
    if previous_good>1
        flag=false;
        for previous_good_index=1:previous_good
            previous_good_field_name=['v' num2str(previous_good_index)];
            if ~isfield(othergoodtest,previous_good_field_name) || ~isfield(othergoodtype,previous_good_field_name), continue; end
            lowerlimit=0.9*othergoodtest.(previous_good_field_name);
            upperlimit=1.1*othergoodtest.(previous_good_field_name);
            if test>lowerlimit && test<upperlimit && flag==0
                think=othergoodtype.(previous_good_field_name);
                good=strcmp(type,think);
                flag=true;
                if good, counter_more_good=counter_more_good+1; end
            end
        end
    end

    %Here are tests based on previous bad identifications
    if counter_bad>1
        for previous_bad_index=1:counter_bad
            previous_bad_field_name=['v' num2str(previous_bad_index)];
            if ~isfield(otherbadtest,previous_bad_field_name), continue; end
            lowerlimit=0.9*otherbadtest.(previous_bad_field_name);
            upperlimit=1.1*otherbadtest.(previous_bad_field_name);
            if test>lowerlimit && test<upperlimit
                good=false;
            end
        end
    end

    %Manual classification if automatic tests fail
    if isempty(good)
        if skip_manual_file_check
            %Retain ambiguous frames when unattended checking is selected.
            good=true;
            manual_override=true;
            counter_good=counter_good+1;
            warning('MEGARA:AmbiguousFileRetained', ...
                [files.list{s} ' could not be classified automatically and ' ...
                'was retained because options.skip_manual_file_check=true.'])
        else
            %Open the Good/Bad interface.
            manual_override=false;
            frame_type=type;
            save('test_line.mat','test_line','frame_type')
            uiwait(untitled2)
            loaded_good=load('good.mat');
            if ~isfield(loaded_good,'good'), error('Manual classification did not create variable good in good.mat.'); end
            good=logical(loaded_good.good);
            if good==1
                previous_good=previous_good+1;
                manual_override=true;
                counter_good=counter_good+1;
                othergoodtest.(['v' num2str(previous_good)])=test;
                othergoodtype.(['v' num2str(previous_good)])=type;
            elseif good==0
                counter_bad=counter_bad+1;
                otherbadtest.(['v' num2str(counter_bad)])=test;
            end
        end
    end
end

%check reverse image by looking at max
[~,maximum_line_index]=max(test_line);
if maximum_line_index<2032 && manual_override==0
    good=false;
end

%compares whether it is labelled right.
if good==0
    badfiles=cat(1,badfiles,{fullfile(pwd,[files.list{s} '.fit'])});
    readout=cat(1,readout,{[files.list{s} ' a  ' type ' is bad']});
    save(['bad_' sanitise_file_label(type) files.list{s} '.mat'],'test_line')
end

%Save accepted profile
if good==1
    save(['good_' sanitise_file_label(type) files.list{s} '.mat'],'test_line')
    readout=cat(1,readout,{[files.list{s} ' a  ' type ' is good']});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function file_label=sanitise_file_label(input_label)

%Create safe file label
if isstring(input_label), input_label=char(input_label); end
file_label=regexprep(input_label,'\s+','_');

end
