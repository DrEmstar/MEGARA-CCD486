function [final_delfns,residual,thefit]=leastsquaresfit_spec_w_broad_delfns(wavelength,spec,delfns,vsini,intrinFWHM,epsx)

%Validate inputs
validate_least_squares_delta_function_inputs(wavelength,spec,delfns,vsini,intrinFWHM,epsx);

%everythings in columns
wavelength=wavelength(:);
spec=spec(:);

%delfns must be a 2-column list of wavelength and equivalent width in mA
delfns=sortrows(delfns,1);
delfnwavs=delfns(:,1);

%Define local line profile
centralwave=mean(wavelength);
step=round(median(diff(wavelength))*1000)/1000;
if step<=0, step=median(diff(wavelength)); end

%we have: vsini, intrinFWHM, epsx
result=make_spectral_line(vsini,intrinFWHM,epsx,centralwave,step);

%make line have equivalent width of 1 and invert
line_profile_wavelength_offset=result(:,1);
line_profile_absorption=1-result(:,2);
line_profile_area=sum(line_profile_absorption)*step;
if ~isfinite(line_profile_area) || line_profile_area<=0, error('The broadened line profile has invalid area.'); end
line_profile_absorption=line_profile_absorption/line_profile_area;

%Set fit parameters
x0=delfns(:,2);
lb=zeros(size(x0));
ub=ones(size(x0))*500;

%do the least-squares fitting
options=optimset('Display','off');
[x,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(@(line_strengths,wavelength_axis) delfn_fitting_local(line_strengths,wavelength_axis,delfnwavs,line_profile_wavelength_offset,line_profile_absorption),x0,wavelength,spec,lb,ub,options);

%Store fitted line list and model spectrum
final_delfns=cat(2,delfnwavs(:),x(:));
thefit=delfn_fitting_local(x,wavelength,delfnwavs,line_profile_wavelength_offset,line_profile_absorption);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_least_squares_delta_function_inputs(wavelength,spec,delfns,vsini,intrinFWHM,epsx)

%Check spectrum inputs
if ~isnumeric(wavelength) || ~isvector(wavelength), error('wavelength must be a numeric vector.'); end
if ~isnumeric(spec) || ~isvector(spec), error('spec must be a numeric vector.'); end
if numel(wavelength)~=numel(spec), error('wavelength and spec must have the same number of elements.'); end
if any(~isfinite(wavelength(:))) || any(~isfinite(spec(:))), error('wavelength and spec must contain only finite values.'); end
if any(diff(wavelength(:))<=0), error('wavelength must be strictly increasing.'); end

%Check delta-function line list
if ~isnumeric(delfns) || size(delfns,2)~=2, error('delfns must be a numeric two-column matrix: wavelength and equivalent width.'); end
if isempty(delfns), error('delfns is empty.'); end
if any(~isfinite(delfns(:))), error('delfns must contain only finite values.'); end

%Check broadening inputs
if ~isnumeric(vsini) || ~isscalar(vsini) || ~isfinite(vsini) || vsini<=0, error('vsini must be a finite scalar value > 0.'); end
if ~isnumeric(intrinFWHM) || ~isscalar(intrinFWHM) || ~isfinite(intrinFWHM) || intrinFWHM<=0, error('intrinFWHM must be a finite scalar value > 0.'); end
if ~isnumeric(epsx) || ~isscalar(epsx) || ~isfinite(epsx) || epsx<0 || epsx>1, error('epsx must be a finite scalar value between 0 and 1.'); end

end

function model_spectrum=delfn_fitting_local(line_strengths,wavelength_axis,delfnwavs,line_profile_wavelength_offset,line_profile_absorption)

%Build fitted spectrum from shifted line profiles
wavelength_axis=wavelength_axis(:);
line_strengths=line_strengths(:);
model_absorption=zeros(size(wavelength_axis));
for line_index=1:numel(delfnwavs)
    shifted_line_wavelength=delfnwavs(line_index)+line_profile_wavelength_offset;
    line_absorption=interp1(shifted_line_wavelength,line_profile_absorption,wavelength_axis,'linear',0);
    model_absorption=model_absorption+(line_strengths(line_index)/1000)*line_absorption;
end
model_spectrum=1-model_absorption;

end

function result=make_spectral_line(vsini,intrinFWHM,epsx,centralwave,step)
%input parameters are:
%vsini        vsini of synth line in km/s
%intrinFWHM   intrinsic FWHM of synth line in km/s
%epsx         limb-darkening parameter
%centralwave  central wavelength used to obtain profile being measured
%outputwave   wavelength vector to output (project) the synth line on to

%Set constants
c=299792.458; %speed of light (km/s)

% Define rotational broadening curve over fine wavelength grid initially
deltalamda=-10:0.001:10;

%prepare parameters for broadening curve
lamda=centralwave;
lam=deltalamda/lamda;
lamgrp=lam/(vsini/c);

%rotational broadening fn from Gray 1992
G=real((2*(1-epsx)*((1-lamgrp.^2).^0.5)+0.5*pi*epsx*(1-lamgrp.^2))/(pi*lamda*vsini/c*(1-epsx/3)));
indices=find(G>0);  %indices to the positive part of G
if isempty(indices), error('Rotational broadening curve has no positive finite region.'); end
if numel(indices)/2==round(numel(indices)/2)  %if even number of elements left in G
    %must remove the one closest to zero as we want an odd number of elements
    ttt=find([G(indices(1)) G(indices(end))]==min([G(indices(1)) G(indices(end))]));
    if ttt==1
        indices(1)=[];
    else
        indices(end)=[];
    end
end
broadcurve=G(min(indices):max(indices));%select the right parts of G

% make gaussian curve to be broadened in same wavelength space as broadening profile
a1=1;
b1=0;

%convert FWHM (velocity to wavelength)
FWHM=intrinFWHM;
FWHM=FWHM*lamda/c;

%make the gaussian
gaussiancurve=a1*exp(-(2*sqrt(log(2))*(deltalamda-b1)/FWHM).^2);

% convolve the gaussian line and the broadening curve to produce the
% broadened curve
res=real(conv(broadcurve,gaussiancurve));

% scale to height 1 and invert
res=res/max(res);
res=1-res;

% make x axis for result
%ceil(numel(res)/2); is the central element
newax=1:numel(res); %create
newax=newax-ceil(numel(res)/2); %centre
newax=newax*0.001; %scale

%trim to just at line limits
d=find(res<0.999);
if isempty(d), error('Could not identify useful limits of broadened line profile.'); end
trim_start=max(1,min(d)-100);     %100 pixels is 0.1A
trim_stop=min(numel(res),max(d)+100);
res=res(trim_start:trim_stop);
newax=newax(trim_start:trim_stop);

%make axis based on step input
st=ceil(newax(1)/step)*step;
ed=floor(newax(end)/step)*step;
finax=st:step:ed;
if isempty(finax), error('Output wavelength axis for broadened line profile is empty.'); end

%project on to requested wavelength axis
newres=spline(newax,res,finax);

%combine xax and int vector
result=finax(:);
result(:,2)=newres(:);

end