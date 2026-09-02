function headers=get_required_headers_from_header(header)

%Initialise required headers
headers.NAXIS1=[];
headers.NAXIS2=[];
headers.HERCEXPT=[];    %Hercules Exposure type
headers.HERCFIB=[];     %Hercules Fibre 1/2/3
headers.DATE_OBS=[];    %Observation date OLD CCD SOFTWARE
headers.DATE=[];        %Observation date NEW CCD SOFTWARE 
headers.MJD_OBS=[];     %Modfd, corr'd, heliocentric, mid obs JD
headers.REC_STRT=[];    %Recorded UT time of observation start
headers.EXPTIME=[];     %Actual exposure time in seconds
headers.HERCFWMT=[];    %Hercules flux weighted mean exp time (mins) ->actually seconds I think ...
headers.HERCFTC=[];     %Hercules flux meter total counts
headers.PORT=[];        %number of ports -1
headers.RA=[];          %RA J2000 coords
headers.DEC=[];         %DEC J2000 coords
headers.RA_2000=[];
headers.DEC_2000=[];
headers.OBSERVER=[];
%following are HERCULES processed headers (using HRSP: prepare_descriptors)
headers.RVCORR=[];      %
headers.JDCORR=[];      %
headers.JD_MID=[];      %
headers.OBJECT=[];      %

%Validate header input
if ~iscell(header), error('header must be a cell array of FITS header cards.'); end

%Read header cards
for header_index=1:numel(header)
    header_card=char(header{header_index});
    header_keyword=get_fits_keyword(header_card);
    header_value=get_fits_value(header_card);
    if strcmp(header_keyword,'NAXIS1')
        headers.NAXIS1=str2double(header_value);
    elseif strcmp(header_keyword,'NAXIS2')
        headers.NAXIS2=str2double(header_value);
    elseif strcmp(header_keyword,'HERCEXPT')
        cleaned_value=clean_fits_string_value(header_value);
        if contains(cleaned_value,'Thorium','IgnoreCase',true)
            headers.HERCEXPT='Thorium';
        elseif contains(cleaned_value,'Stellar','IgnoreCase',true)
            headers.HERCEXPT='Stellar';
        elseif contains(cleaned_value,'White','IgnoreCase',true) || contains(cleaned_value,'Smooth','IgnoreCase',true)
            headers.HERCEXPT='White L';
        end
    elseif strcmp(header_keyword,'HERCFIB')
        headers.HERCFIB=parse_hercules_fibre(header_card,header_value);
    elseif strcmp(header_keyword,'DATE-OBS')
        headers.DATE_OBS=clean_fits_string_value(header_value);
    elseif strcmp(header_keyword,'DATE')         % NEW CCD SOFTWARE
        headers.DATE=clean_fits_string_value(header_value);                                % NEW CCD SOFTWARE
    elseif strcmp(header_keyword,'MJD-OBS')
        headers.MJD_OBS=str2double(header_value);
    elseif strcmp(header_keyword,'REC-STRT')
        headers.REC_STRT=clean_fits_string_value(header_value);
    elseif strcmp(header_keyword,'EXPTIME')
        headers.EXPTIME=str2double(header_value);
    elseif strcmp(header_keyword,'HERCFWMT')
        headers.HERCFWMT=str2double(header_value);
    elseif strcmp(header_keyword,'HERCFTC')
        headers.HERCFTC=str2double(header_value);
%   elseif ~isempty(strfind(header{m},'PARAM59')) 
    elseif strcmp(header_keyword,'PARAM58') % for NEW CCD SOFTWARE   
        headers.PORT=parse_port_number(header,header_index,header_value); % for OLD CCD SOFTWARE
    elseif strcmp(header_keyword,'RVCORR')
        headers.RVCORR=str2double(header_value);
    elseif strcmp(header_keyword,'JDCORR')
        headers.JDCORR=str2double(header_value);
    elseif strcmp(header_keyword,'JD-MID')
        headers.JD_MID=str2double(header_value);
    elseif strcmp(header_keyword,'OBSERVER')
        headers.OBSERVER=clean_fits_string_value(header_value);
    elseif contains(header_card,'COMMENT  RA')
        headers.RA_2000=parse_comment_coordinate(header_card,15,24);
    elseif contains(header_card,'COMMENT  DEC')
        headers.DEC_2000=parse_comment_coordinate(header_card,16,24);
    elseif strcmp(header_keyword,'OBJECT')
        headers.OBJECT=clean_fits_string_value(header_value);
    end
end

%check for important missing headers!
if isempty(headers.MJD_OBS)
    if isempty(headers.JD_MID)
        if isempty(headers.DATE_OBS) && isempty(headers.EXPTIME) && isempty(headers.REC_STRT)
            warning('NOTE: there is insufficient header information to determine JD of observation!!!!!!!')
        end
    end
end

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

function fibre_number=parse_hercules_fibre(header_card,header_value)

%Parse fibre number
fibre_number=str2double(header_value);
if isnan(fibre_number)
    if numel(header_card)>=12 && header_card(12)~=' '
        fibre_number=str2double(header_card(12));
    elseif numel(header_card)>=28
        fibre_number=str2double(header_card(28));
    else
        fibre_number=[];
    end
end

end

function port_number=parse_port_number(header,header_index,header_value)

%Parse port number
port_value=str2double(header_value);
if isnan(port_value)
    port_number=[];
    return
end
if port_value>10% for OLD CCD SOFTWARE
    if header_index<numel(header)
        next_header_value=get_fits_value(char(header{header_index+1}));
        next_port_value=str2double(next_header_value);
        if isnan(next_port_value)
            port_number=[];
        else
            port_number=next_port_value+1;
        end
    else
        port_number=[];
    end
else
    port_number=port_value+1;
end

end

function coordinate_value=parse_comment_coordinate(header_card,start_index,stop_index)

%Parse comment coordinate
if numel(header_card)>=stop_index
    coordinate_value=strtrim(header_card(start_index:stop_index));
else
    coordinate_value=[];
end

end