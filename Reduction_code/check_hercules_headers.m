function header=check_hercules_headers(head)

%Copy input header
header=head;

%Return if MJD already exists
if isfield(head,'MJD_OBS') && ~isempty(head.MJD_OBS)
    return
end

%Check required timing fields
if ~isfield(head,'DATE_OBS') || isempty(head.DATE_OBS) || ~isfield(head,'REC_STRT') || isempty(head.REC_STRT) || ~isfield(head,'HERCFWMT') || isempty(head.HERCFWMT)
    warning('not enough information to reconstruct MJD time from DATE_OBS REC_STRT and HERCFWMT, exiting ...')
    return
end

%Parse observation date
date=head.DATE_OBS;
date=clean_header_string(date);
if numel(date)<10
    warning('DATE_OBS is too short to reconstruct MJD_OBS.')
    return
end
year=str2double(date(1:4));
month=str2double(date(6:7));
day=str2double(date(9:10));
if any(~isfinite([year,month,day]))
    warning('DATE_OBS could not be parsed to reconstruct MJD_OBS.')
    return
end

%Parse observation start time
time=head.REC_STRT;
time=clean_header_string(time);
if isempty(time)
    warning('REC_STRT is empty and MJD_OBS could not be reconstructed.')
    return
end
if contains(time,'.')
    time=[time '0'];
end

%Calculate start JD
JDstart=convert_date_to_JD(year,month,day,time);

%flux weighted exposure time added
expt=head.HERCFWMT/86400;
if ~isfinite(expt)
    warning('HERCFWMT is not finite and MJD_OBS could not be reconstructed.')
    return
end

%Calculate midpoint MJD
JDmid=JDstart+expt;

%to keep inline with normally written fits headers
JDmid=JDmid-2.4e6;
header.MJD_OBS=JDmid;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output=clean_header_string(input_value)

%Clean quoted header string
if isstring(input_value), input_value=char(input_value); end
if ~ischar(input_value), output=''; return; end
output=strtrim(input_value);
output=strrep(output,'''','');
output=strrep(output,'"','');
output=strtrim(output);

end