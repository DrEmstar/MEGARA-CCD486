function [correction,details]=calculate_barycentric_correction(expdate,midtime_tt,target_or_ra,declination)
%CALCULATE_BARYCENTRIC_CORRECTION HERCULES barycentric RV correction.
%   CORRECTION is in km/s using the existing MEGARA sign convention.
%   EXPDATE is yyyy-MM-dd (a longer ISO date is also
%   accepted), and MIDTIME_TT is HH:mm:ss.SSS in Terrestrial Time.
%
%   CALCULATE_BARYCENTRIC_CORRECTION(DATE,TIME,TARGET) looks up an
%   HD, HIP, HR, Bayer, or Flamsteed target in barycentric_stellar_catalogue.mat and applies the same
%   Hipparcos parallax, proper motion, radial velocity, and 1991.25 epoch
%   conventions as HRSP. If TARGET is not identified, a warning is issued
%   and CORRECTION is NaN.
%
%   CALCULATE_BARYCENTRIC_CORRECTION_MATLAB(DATE,TIME,RA,DEC) retains the
%   original coordinate-only mode, with zero parallax and space motion.
%
%   This implementation reads the DE405 ephemeris distributed with HRSP,
%   but does not use compiled HRSP code or its binary stellar catalogues.

validate_text_input(expdate,'expdate');
validate_text_input(midtime_tt,'midtime_tt');
validate_text_input(target_or_ra,'target_or_ra');
if nargin>=4, validate_text_input(declination,'declination'); end

[year,month,day]=parse_iso_date(expdate);
[hour,minute,second]=parse_clock_time(midtime_tt);
jd_tt=calendar_to_julian_date(year,month,day,hour,minute,second);

%DE405 is evaluated in TDB. This is the same short-period approximation
%used by HRSP.
days_from_j2000=jd_tt-2451545.0;
tdb_minus_tt_seconds=0.001657*sind(357.53+0.98560028*days_from_j2000)+ ...
    0.000022*sind(246.11+0.90251792*days_from_j2000);
jd_tdb=jd_tt+tdb_minus_tt_seconds/86400;

[earth_barycentric_position,earth_barycentric_velocity]= ...
    get_de405_earth_barycentric_state(jd_tdb);

if nargin<4 || isempty(declination)
    [star,identified]=find_barycentric_catalogue_target(target_or_ra);
    if ~identified
        correction=NaN;
        details=struct('target',string(target_or_ra),'identified',false);
        return
    end
    catalogue_mode=true;
else
    star=struct('ra_deg',parse_sexagesimal_coordinate(target_or_ra,true), ...
        'dec_deg',parse_sexagesimal_coordinate(declination,false), ...
        'parallax_mas',0,'pmra_mas_yr',0,'pmdec_mas_yr',0, ...
        'rv_km_s',0,'epoch',2000.0);
    catalogue_mode=false;
end

%Use HRSP's full IAU 2006/2000A Earth orientation and its bundled Delta-T
%table to transform the Earth-fixed HERCULES state into GCRS.
[site_position,site_velocity]=hercules_observatory_state(jd_tt);

observer_barycentric_position=earth_barycentric_position+ ...
    site_position/149597870.691;
source_direction=propagated_direction(star,jd_tdb, ...
    observer_barycentric_position);
orbital_correction=dot(earth_barycentric_velocity,source_direction);

rotational_correction=dot(site_velocity,source_direction);

correction=orbital_correction+rotational_correction;
if ~isfinite(correction), error('Calculated barycentric correction is not finite.'); end

details=struct('identified',true,'catalogue_mode',catalogue_mode, ...
    'source_direction',source_direction,'orbital_correction',orbital_correction, ...
    'rotational_correction',rotational_correction,'star',star);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DE405 interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [earth_position,earth_velocity]=get_de405_earth_barycentric_state(jd_tdb)

persistent ephemeris
if isempty(ephemeris)
    ephemeris_path=fullfile(fileparts(mfilename('fullpath')),'de405');
    ephemeris=read_de405_header(ephemeris_path);
end

if jd_tdb<ephemeris.first_epoch || jd_tdb>ephemeris.last_epoch
    error('Observation epoch %.9f is outside the bundled DE405 range %.1f to %.1f.', ...
        jd_tdb,ephemeris.first_epoch,ephemeris.last_epoch)
end

record_number=floor((jd_tdb-ephemeris.first_epoch)/ ...
    ephemeris.days_per_record)+2;
if jd_tdb==ephemeris.last_epoch, record_number=record_number-1; end
record_start=(record_number-2)*ephemeris.days_per_record+ ...
    ephemeris.first_epoch;
normalised_time=(jd_tdb-record_start)/ephemeris.days_per_record;

file_identifier=fopen(ephemeris.path,'r','ieee-le');
if file_identifier==-1, error('Could not open DE405 ephemeris: %s',ephemeris.path); end
cleanup_file=onCleanup(@() fclose(file_identifier));
fseek(file_identifier,record_number*ephemeris.record_size,'bof');
record=fread(file_identifier,ephemeris.record_size/8,'double=>double');
if numel(record)~=ephemeris.record_size/8
    error('Could not read the requested DE405 record for JD %.9f.',jd_tdb)
end

%C/DE405 target indices are EMBARY=2, MOON=9, SUN=10. The Sun is not
%needed for barycentric Earth velocity, but the numbering is retained here
%for direct comparison with the HRSP implementation.
[embary_position,embary_velocity]=interpolate_de405_target(record,normalised_time, ...
    ephemeris,2);
[moon_position,moon_velocity]=interpolate_de405_target(record,normalised_time, ...
    ephemeris,9);

earth_position=embary_position-moon_position/(1+ephemeris.earth_moon_mass_ratio);
earth_velocity_au_per_day=embary_velocity- ...
    moon_velocity/(1+ephemeris.earth_moon_mass_ratio);
au_per_day_to_km_per_second=1731.45683670138888888889;
earth_velocity=earth_velocity_au_per_day*au_per_day_to_km_per_second;

end


function header=read_de405_header(ephemeris_path)

if ~isfile(ephemeris_path), error('Bundled DE405 ephemeris not found: %s',ephemeris_path); end
file_identifier=fopen(ephemeris_path,'r','ieee-le');
if file_identifier==-1, error('Could not open DE405 ephemeris: %s',ephemeris_path); end
cleanup_file=onCleanup(@() fclose(file_identifier));

header=struct();
header.path=ephemeris_path;
header.first_epoch=read_binary_scalar(file_identifier,2652,'double');
header.last_epoch=read_binary_scalar(file_identifier,2660,'double');
header.days_per_record=read_binary_scalar(file_identifier,2668,'double');
header.astronomical_unit=read_binary_scalar(file_identifier,2680,'double');
header.earth_moon_mass_ratio=read_binary_scalar(file_identifier,2688,'double');
header.pointer=zeros(13,1);
header.coefficient_count=zeros(13,1);
header.subinterval_count=zeros(13,1);

for target_index=0:11
    offset=2696+12*target_index;
    header.pointer(target_index+1)=read_binary_scalar(file_identifier,offset,'int32');
    header.coefficient_count(target_index+1)=read_binary_scalar(file_identifier,offset+4,'int32');
    header.subinterval_count(target_index+1)=read_binary_scalar(file_identifier,offset+8,'int32');
end
header.pointer(13)=read_binary_scalar(file_identifier,2844,'int32');
header.coefficient_count(13)=read_binary_scalar(file_identifier,2848,'int32');
header.subinterval_count(13)=read_binary_scalar(file_identifier,2852,'int32');

[maximum_pointer,maximum_index]=max(header.pointer);
if maximum_index==12
    number_dimensions=2;
else
    number_dimensions=3;
end
kernel_size=2*(maximum_pointer+number_dimensions* ...
    header.coefficient_count(maximum_index)* ...
    header.subinterval_count(maximum_index)-1);
header.record_size=kernel_size*4;

end


function value=read_binary_scalar(file_identifier,offset,precision)

if fseek(file_identifier,offset,'bof')~=0, error('Could not seek within DE405 ephemeris.'); end
value=fread(file_identifier,1,[precision '=>double']);
if isempty(value), error('Could not read DE405 ephemeris header.'); end

end


function [position,velocity]=interpolate_de405_target(record,normalised_time,header,target_index)

array_index=target_index+1;
number_coefficients=header.coefficient_count(array_index);
number_subintervals=header.subinterval_count(array_index);
scaled_time=number_subintervals*normalised_time;
subinterval=floor(scaled_time);
if subinterval>=number_subintervals, subinterval=number_subintervals-1; end
chebyshev_time=2*(scaled_time-subinterval)-1;

position_polynomial=zeros(number_coefficients,1);
velocity_polynomial=zeros(number_coefficients,1);
position_polynomial(1)=1;
if number_coefficients>=2
    position_polynomial(2)=chebyshev_time;
    velocity_polynomial(2)=1;
end
if number_coefficients>=3, velocity_polynomial(3)=4*chebyshev_time; end
for coefficient_index=3:number_coefficients
    position_polynomial(coefficient_index)=2*chebyshev_time* ...
        position_polynomial(coefficient_index-1)- ...
        position_polynomial(coefficient_index-2);
end
for coefficient_index=4:number_coefficients
    velocity_polynomial(coefficient_index)=2*chebyshev_time* ...
        velocity_polynomial(coefficient_index-1)+ ...
        2*position_polynomial(coefficient_index-1)- ...
        velocity_polynomial(coefficient_index-2);
end

position=zeros(1,3);
velocity=zeros(1,3);
target_start=header.pointer(array_index);
subinterval_size=3*number_coefficients;
velocity_scale=2*number_subintervals/header.days_per_record;
for component_index=1:3
    coefficient_start=target_start+subinterval*subinterval_size+ ...
        (component_index-1)*number_coefficients;
    coefficients=record(coefficient_start: ...
        coefficient_start+number_coefficients-1);
    position(component_index)=dot(position_polynomial,coefficients)/ ...
        header.astronomical_unit;
    velocity(component_index)=dot(velocity_polynomial,coefficients)* ...
        velocity_scale/header.astronomical_unit;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coordinates, time, and observatory motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function direction=propagated_direction(star,jd_tdb,observer_position_au)

ra=deg2rad(star.ra_deg);
dec=deg2rad(star.dec_deg);
p=[cos(dec)*cos(ra),cos(dec)*sin(ra),sin(dec)];
e_ra=[-sin(ra),cos(ra),0];
e_dec=[-sin(dec)*cos(ra),-sin(dec)*sin(ra),cos(dec)];
dt_years=(jd_tdb-(2451545.0+(star.epoch-2000.0)*365.25))/365.25;

parallax_arcsec=star.parallax_mas*1e-3;
pmra_arcsec_yr=star.pmra_mas_yr*1e-3;
pmdec_arcsec_yr=star.pmdec_mas_yr*1e-3;

if parallax_arcsec>0
    distance_pc=1/parallax_arcsec;
    transverse_constant=4.74047046324815575328922;
    velocity=(transverse_constant*pmra_arcsec_yr/parallax_arcsec)*e_ra+ ...
        (transverse_constant*pmdec_arcsec_yr/parallax_arcsec)*e_dec+ ...
        star.rv_km_s*p;
    year_seconds=31557600.0;
    parsec_km=30856775813057.2895329155;
    position_pc=distance_pc*p+velocity*(dt_years*year_seconds/parsec_km);
    parsec_au=206264.80624709636;
    position_pc=position_pc-observer_position_au/parsec_au;
    direction=position_pc/norm(position_pc);
else
    proper_motion=pmra_arcsec_yr*e_ra+pmdec_arcsec_yr*e_dec;
    rate=norm(proper_motion);
    if rate==0
        direction=p;
    else
        axis=cross(p,proper_motion);
        axis=axis/norm(axis);
        angle=deg2rad(rate*dt_years/3600);
        %HRSP's arbitrary_rotation rotates the reference frame, so its
        %positive angle has the opposite sign to active Rodrigues rotation.
        direction=p*cos(angle)-cross(axis,p)*sin(angle)+ ...
            axis*dot(axis,p)*(1-cos(angle));
    end
    %Retain HRSP's cent_obs behaviour for its proper-motion-only model.
    direction=direction-observer_position_au/206264.80624709636;
    direction=direction/norm(direction);
end
end


function degrees=parse_sexagesimal_coordinate(value,is_right_ascension)

value=clean_text_value(value);
negative=startsWith(strtrim(value),'-');
value=strrep(value,':',' ');
parts=sscanf(value,'%f');
if isempty(parts), error('Could not parse coordinate: %s',value); end

if numel(parts)==1
    degrees=parts(1);
    if is_right_ascension && abs(degrees)<=24, degrees=15*degrees; end
    return
end
if numel(parts)<3, parts(3)=0; end
magnitude=abs(parts(1))+abs(parts(2))/60+abs(parts(3))/3600;
if negative || parts(1)<0, magnitude=-magnitude; end
if is_right_ascension, magnitude=15*magnitude; end
degrees=magnitude;

end


function jd=calendar_to_julian_date(year,month,day,hour,minute,second)

if month<=2, year=year-1; month=month+12; end
century=floor(year/100);
calendar_term=2-century+floor(century/4);
jd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+day+ ...
    calendar_term-1524.5+(hour+minute/60+second/3600)/24;

end


function [year,month,day]=parse_iso_date(value)

value=clean_text_value(value);
tokens=regexp(value,'^(\d{4})-(\d{2})-(\d{2})','tokens','once');
if isempty(tokens), error('expdate must start with yyyy-MM-dd.'); end
year=str2double(tokens{1});
month=str2double(tokens{2});
day=str2double(tokens{3});

end


function [hour,minute,second]=parse_clock_time(value)

value=clean_text_value(value);
tokens=regexp(value,'(\d{1,2}):(\d{2}):([0-9.]+)','tokens','once');
if isempty(tokens), error('midtime_tt must contain HH:mm:ss.SSS.'); end
hour=str2double(tokens{1});
minute=str2double(tokens{2});
second=str2double(tokens{3});

end


function validate_text_input(value,label)

if ~(ischar(value) || (isstring(value) && isscalar(value)))
    error('%s must be a character vector or string scalar.',label)
end
if strlength(string(value))==0, error('%s must not be empty.',label); end

end


function value=clean_text_value(value)

value=strtrim(char(value));
if numel(value)>=2 && value(1)=='"' && value(end)=='"'
    value=value(2:end-1);
end
if numel(value)>=2 && value(1)=='''' && value(end)==''''
    value=value(2:end-1);
end
value=strtrim(value);

end
