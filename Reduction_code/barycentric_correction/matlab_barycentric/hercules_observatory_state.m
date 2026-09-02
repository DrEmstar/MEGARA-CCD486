function [position_gcrs,velocity_gcrs,earth_to_gcrs]=hercules_observatory_state(jd_tt)
%HERCULES_OBSERVATORY_STATE Calculate HERCULES topocentric state vectors.
%   Position is in km and velocity in km/s, both in GCRS/ICRS axes.

persistent data
if isempty(data)
    table_path=fullfile(fileparts(mfilename('fullpath')), ...
        'barycentric_earth_orientation.mat');
    if ~isfile(table_path)
        error(['HRSP Earth-orientation table not found: %s\nRun ' ...
            'convert_hrsp_earth_orientation_to_mat first.'],table_path)
    end
    loaded=load(table_path,'orientation');
    data=loaded.orientation;
end

tt_minus_utc=terrestrial_time_minus_utc(jd_tt);
jd_utc=jd_tt-tt_minus_utc/86400;
deltat=interpolate_hrsp_deltat(jd_utc,data);
jd_ut1=jd_tt-deltat/86400;
earth_to_gcrs=earth_orientation_matrix(jd_tt,jd_ut1,data);

longitude=deg2rad(15*(11+21/60+51.6/3600));
latitude=deg2rad(-(43+59/60+12/3600));
height_km=1.027;
equatorial_radius=6378.137;
flattening=1/298.257223563;
q=(1-flattening)^2;
C=1/sqrt(cos(latitude)^2+q*sin(latitude)^2);
S=q*C;
position_ecef=[(equatorial_radius*C+height_km)*cos(latitude)*cos(longitude); ...
    (equatorial_radius*C+height_km)*cos(latitude)*sin(longitude); ...
    (equatorial_radius*S+height_km)*sin(latitude)];

era_coeff=0.00273781191135448;
angular_velocity=(1+era_coeff)*2*pi/86400;
velocity_ecef=cross([0;0;angular_velocity],position_ecef);
position_gcrs=(earth_to_gcrs*position_ecef).';
velocity_gcrs=(earth_to_gcrs*velocity_ecef).';
end


function matrix=earth_orientation_matrix(jd_tt,jd_ut1,data)

t=(jd_tt-2451545.0)/36525;
fw=zeros(1,4);
for row=1:4, fw(row)=arcsec_to_rad(polyval(fliplr(data.fukushima_williams(row,:)),t)); end
gcrs_to_true=fw_matrix(fw);
[dpsi,deps]=nut06a(t,data);
fw(3)=fw(3)+dpsi;
fw(4)=fw(4)+deps;
gcrs_to_true=fw_matrix(fw);

x=gcrs_to_true(3,1);
y=gcrs_to_true(3,2);
z=gcrs_to_true(3,3);
rr=x*x+y*y;
if rr>0, e=atan2(y,x); else, e=0; end
d=atan(sqrt(rr/(1-rr)));
s=cio_locator(t,x,y,data);

a=1/(1+z);
rs=[1-a*x*x,-a*x*y,-x];
p=dot(gcrs_to_true(1,:),rs);
q=dot(gcrs_to_true(2,:),rs);
equation_origin=s;
if p~=0 || q~=0, equation_origin=equation_origin-atan2(q,p); end

djd=jd_ut1-2451545.0;
jd0=floor(jd_ut1+0.5)-0.5;
jdf=jd_ut1-jd0;
era_fraction=0.7790572732640+rem(0.00273781191135448*djd,1)+jdf+0.5;
earth_rotation_angle=mod(era_fraction,1)*2*pi;
gast=earth_rotation_angle-equation_origin;
gcrs_to_earth=rotation_matrix(3,gast)*gcrs_to_true;
matrix=gcrs_to_earth.';

%Keep the HRSP construction above explicit; e and d are also the angles
%used for its CIO matrix, which is not otherwise needed by set_topo.
if ~isfinite(e+d), error('Invalid HRSP celestial-pole orientation.'); end
end


function matrix=fw_matrix(fw)

matrix=rotation_matrix(1,-fw(4))*rotation_matrix(3,-fw(3))* ...
    rotation_matrix(1,fw(2))*rotation_matrix(3,fw(1));
end


function [dpsi,deps]=nut06a(t,data)

fundamental=zeros(5,1);
for j=1:5
    fundamental(j)=rem(polyval(fliplr(data.fundamental_arguments(j,:)),t), ...
        1296000)*4.848136811095359935899141e-6;
end
arguments=data.nals*fundamental;
sine=sin(arguments);
cosine=cos(arguments);
dpsils=sum((data.cls(:,1)+data.cls(:,2)*t).*sine+data.cls(:,3).*cosine);
depsls=sum((data.cls(:,4)+data.cls(:,5)*t).*cosine+data.cls(:,6).*sine);

planetary=zeros(14,1);
for j=1:14
    value=data.planetary_arguments(j,1)+data.planetary_arguments(j,2)*t;
    if j<14, planetary(j)=rem(value,2*pi); else, planetary(j)=value*t; end
end
arguments=data.napl*planetary;
sine=sin(arguments);
cosine=cos(arguments);
dpsipl=sum(data.cpl(:,1).*sine+data.cpl(:,2).*cosine);
depspl=sum(data.cpl(:,3).*sine+data.cpl(:,4).*cosine);

scale=4.848136811095359935899141e-13;
dpsi=(dpsils+dpsipl)*scale;
deps=(depsls+depspl)*scale;
factor=-2.7774e-6*t;
dpsi=dpsi+dpsi*(factor+0.4697e-6);
deps=deps+deps*factor;
end


function s=cio_locator(t,x,y,data)

fa=zeros(8,1);
for j=1:5
    fa(j)=rem(polyval(fliplr(data.fundamental_arguments(j,:)),t), ...
        1296000)*4.848136811095359935899141e-6;
end
fa(6)=rem(data.planetary_arguments(6,1)+data.planetary_arguments(6,2)*t,2*pi);
fa(7)=rem(data.planetary_arguments(7,1)+data.planetary_arguments(7,2)*t,2*pi);
value=data.planetary_arguments(14,1)+data.planetary_arguments(14,2)*t;
fa(8)=value*t;

coefficients=data.cio_sp;
for order=1:5
    first=data.cio_ia(order)+1;
    last=data.cio_ib(order)+1;
    arguments=data.cio_ks(first:last,:)*fa;
    coefficients(order)=coefficients(order)+sum( ...
        data.cio_ss(first:last,1).*sin(arguments)+ ...
        data.cio_ss(first:last,2).*cos(arguments));
end
s=polyval(fliplr(coefficients),t);
s=arcsec_to_rad(s)-x*y/2;
end


function matrix=rotation_matrix(axis,angle)

c=cos(angle); s=sin(angle);
switch axis
    case 1, matrix=[1 0 0;0 c s;0 -s c];
    case 2, matrix=[c 0 -s;0 1 0;s 0 c];
    case 3, matrix=[c s 0;-s c 0;0 0 1];
    otherwise, error('Invalid rotation axis.');
end
end


function seconds=interpolate_hrsp_deltat(jd_utc,data)

if jd_utc<data.deltat_jd(1)
    seconds=data.deltat_seconds(1);
elseif jd_utc>=data.deltat_jd(end)
    seconds=data.deltat_seconds(end);
else
    lower=find(data.deltat_jd<=jd_utc,1,'last');
    upper=lower+1;
    fraction=(jd_utc-data.deltat_jd(lower))/ ...
        (data.deltat_jd(upper)-data.deltat_jd(lower));
    seconds=data.deltat_seconds(lower)+fraction* ...
        (data.deltat_seconds(upper)-data.deltat_seconds(lower));
end
end


function seconds=terrestrial_time_minus_utc(jd_tt)

effective_jd=[2441317.5 2441499.5 2441683.5 2442048.5 2442413.5 ...
    2442778.5 2443144.5 2443509.5 2443874.5 2444239.5 2444786.5 ...
    2445151.5 2445516.5 2446247.5 2447161.5 2447892.5 2448257.5 ...
    2448804.5 2449169.5 2449534.5 2450083.5 2450630.5 2451179.5 ...
    2453736.5 2454832.5 2456109.5 2457204.5 2457754.5];
tai_minus_utc=10:37;
utc_estimate=jd_tt-(tai_minus_utc(end)+32.184)/86400;
active=find(effective_jd<=utc_estimate,1,'last');
if isempty(active), error('Earth orientation currently supports dates from 1972 onward.'); end
seconds=tai_minus_utc(active)+32.184;
end


function radians=arcsec_to_rad(arcseconds)

radians=arcseconds*4.848136811095359935899141e-6;
end
