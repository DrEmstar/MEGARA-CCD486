function summed_flat=merge_flats(info,Bias)

%find all flat fields and check that they do not saturate and that they all agree with each other and then add them

%Do you want to ignore flats that are within 5000 ADU of saturation?')
ignore_near_saturated_flats=false;%input('(y/[n]) ->','s');

%Set saturation threshold
maxdatalimit=65526;
near_saturation_margin=5000;

%Initialise output
summed_flat=[];
number_of_flats_used=0;

%Check file list
if ~isstruct(info) || ~isfield(info,'list'), error('info must be a structure containing info.list.'); end
numfiles=numel(info.list);

%Sum valid flat fields
for file_index=1:numfiles
    file_identifier=info.list{file_index};
    if ~isfield(info,file_identifier) || ~isfield(info.(file_identifier),'HERCEXPT'), continue; end
    temptype=info.(file_identifier).HERCEXPT;
    if file_index==1
        fprintf('Reading and merging flat frames...\n')
    end
    if strcmp(temptype,'White L')
        %fprintf('reading file %s.fit\n',file_identifier)
        try
            data=raw_hercules_fitsread([file_identifier '.fit']);
        catch read_error
            warning(['Could not read flat file ' file_identifier '.fit: ' read_error.message])
            continue
        end
        if ignore_near_saturated_flats && max(data(:))>maxdatalimit-near_saturation_margin
            continue
        else
            if ~isempty(Bias)
                if isequal(size(Bias),size(data))
                    data=data-Bias;
                else
                    warning(['Bias frame size does not match flat file ' file_identifier '.fit. Bias subtraction skipped for this flat.'])
                end
            end
            number_of_flats_used=number_of_flats_used+1;
            if number_of_flats_used==1
                summed_flat=data;
            else
                summed_flat=summed_flat+data;
            end
        end
    end
end

%Warn if no flats were available
if number_of_flats_used==0
    warning('No valid White L flat fields were found for flat merging.')
    summed_flat=[];
end

%Report completion
disp('flat summing complete')

end