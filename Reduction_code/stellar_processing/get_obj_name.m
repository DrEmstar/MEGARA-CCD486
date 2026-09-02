function obj_name=get_obj_name(filename)

%Initialise output
obj_name='';

%Read FITS header
header=get_header(filename);

%Search for OBJECT header
for header_index=1:numel(header)
    header_card=char(header{header_index});
    header_keyword=get_fits_keyword(header_card);
    if strcmp(header_keyword,'OBJECT')
        obj_name=get_fits_value(header_card);
        obj_name=clean_fits_string_value(obj_name);
        return
    end
end

%Warn if OBJECT is missing
disp('no OBJECT header found!')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function header_keyword=get_fits_keyword(header_card)

%Extract FITS keyword
if numel(header_card)<8
    header_keyword=strtrim(header_card);
else
    header_keyword=strtrim(header_card(1:8));
end

end

function header_value=get_fits_value(header_card)

%Extract FITS value
header_value='';
equals_index=find(header_card=='=',1,'first');
if isempty(equals_index), return; end
value_start=equals_index+1;
comment_index=find(header_card(value_start:end)=='/',1,'first');
if isempty(comment_index)
    value_stop=numel(header_card);
else
    value_stop=value_start+comment_index-2;
end
header_value=strtrim(header_card(value_start:value_stop));

end

function cleaned_value=clean_fits_string_value(header_value)

%Clean FITS string value
cleaned_value=strtrim(header_value);
cleaned_value=strrep(cleaned_value,'''','');
cleaned_value=strtrim(cleaned_value);

end