function JD=convert_date_to_JD(year,month,day,time)
%only for obs taken after the start of the year 2000!
%time in ??:??:??.?? format entered as a string - MUST be UT 24 hour time
%year month day entered as numbers

%The Julian date for CE  2000 January  1st  00:00:00.0 UT is
%JD 2451544.50000

%Validate date inputs
validate_date_inputs(year,month,day,time);

%Normalise year
if year<100 %two digit year was entered
    year=year+2000;
end

%Parse time string
[hours,minutes,seconds]=parse_time_string(time);

%Calculate day fraction
totalsecs=hours*3600+minutes*60+seconds;
dayfraction=totalsecs/86400; % 86400 seconds in a day close approx.

%Convert calendar date to Julian Date
if month<=2
    year=year-1;
    month=month+12;
end
century=floor(year/100);
gregorian_correction=2-century+floor(century/4);
JD=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+gregorian_correction-1524.5+dayfraction;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_date_inputs(year,month,day,time)

%Check date values
if ~isnumeric(year) || ~isscalar(year) || ~isfinite(year), error('year must be a finite scalar numeric value.'); end
if ~isnumeric(month) || ~isscalar(month) || ~isfinite(month) || month<1 || month>12, error('month must be a finite scalar numeric value from 1 to 12.'); end
if ~isnumeric(day) || ~isscalar(day) || ~isfinite(day) || day<1 || day>31, error('day must be a finite scalar numeric value from 1 to 31.'); end
if isstring(time), time=char(time); end
if ~ischar(time) || isempty(time), error('time must be a character vector or string scalar.'); end

end

function [hours,minutes,seconds]=parse_time_string(time)

%Clean time string
if isstring(time), time=char(time); end
time=strtrim(time);
time=strrep(time,'''','');
time=strrep(time,'"','');
time_parts=strsplit(time,':');

if numel(time_parts)<3
    error('time must be in HH:MM:SS or HH:MM:SS.SSS format.')
end

%Parse time fields
hours=str2double(time_parts{1});
minutes=str2double(time_parts{2});
seconds=str2double(time_parts{3});

if any(~isfinite([hours,minutes,seconds]))
    error('time contains invalid numeric values.')
end

if hours<0 || hours>=24
    error('hours must be from 0 to less than 24.')
end

if minutes<0 || minutes>=60
    error('minutes must be from 0 to less than 60.')
end

%Allow seconds=60 from rounded headers, then normalise
if seconds<0 || seconds>60
    error('seconds must be from 0 to 60.')
end

if seconds==60
    seconds=0;
    minutes=minutes+1;
end

if minutes==60
    minutes=0;
    hours=hours+1;
end

%Allow 23:59:60 to become 24:00:00.
%The main JD calculation handles this because dayfraction becomes 1.
if hours==24
    hours=0;
    seconds=seconds+86400;
end

end