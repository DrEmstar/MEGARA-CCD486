function [vsini,fittedint,int]=vsini_directfit(vel,int,locs)
%format long g
%global lamda0 epsx
% Fit the line profile by optimizing an unconstrained FWHM Gaussian convolved with the rotational broadening curve of Gray 1992 for the specified wavelength
% (usually the weighted centre wavelength of the ccf-function). Note that the input spectrum is in velocity space.

%Validate inputs
validate_vsini_directfit_inputs(vel,int,locs);

%Preserve input orientation
return_column=iscolumn(vel);
vel=vel(:)';
int=int(:)';

%Initialise constants
c=299792.458; % Speed of light

%Select line centre
line_centre_velocity=locs(1);

%invert line profile and normalise
int=1-int;
maximum_line_depth=max(int);
if ~isfinite(maximum_line_depth) || maximum_line_depth<=0
    vsini=NaN;
    fittedint=NaN(size(int));
    int=1-int;
    if return_column, fittedint=fittedint(:); int=int(:); end
    return
end
int=int/maximum_line_depth;

%function parameters
a1=line_centre_velocity;%input line centre (in km/s), this will be subtracted from the velocities to centre the line
vel=vel-a1;
a3=10;%input starting V sin(i) (in km/s, a range will be tested)
a2=a3/2;%input starting FWHM for pre-rotationally broadened line (in km/s, a range will be tested)
lambda0=5500;%input lambda0 (use centre wavelength of ccf if a ccf)
limb_darkening_coefficient=0.555;%input limb darkening coefficient

%Set fit bounds
lower_bounds=[0 0];
upper_bounds=[200 1000];

%respfnwidth=input('input spectrograph response width parameter or leave empty to fit Thar lines (telwave and telint); -> ');

%convert FWHM to wavelength
initial_parameters=[a2,a3];
initial_parameters(1)=a2(1)*lambda0/c;

%fit line data with a broadened Gaussian
fit_options=optimset('Display','off');
try
    [fit_parameters,resnorm,residual,exitflag]=lsqcurvefit(@(fit_parameters,velocity_axis) vsini_fittingfn(fit_parameters,velocity_axis,lambda0,limb_darkening_coefficient),initial_parameters,vel,int,lower_bounds,upper_bounds,fit_options);
catch
    vsini=NaN;
    fittedint=NaN(size(int));
    int=1-int;
    if return_column, fittedint=fittedint(:); int=int(:); end
    return
end

%Convert fitted width back to velocity
results=fit_parameters;
fit_parameters(1)=fit_parameters(1)/lambda0*c;
results1=fit_parameters;

%Extract fitted values
vsini=results1(2);
fittedint=vsini_fittingfn(results,vel,lambda0,limb_darkening_coefficient);
fittedint=1-fittedint;
int=1-int;

%Restore output orientation
if return_column
    fittedint=fittedint(:);
    int=int(:);
end

% figure
% plot(vel,1-int,'.k')
% hold on
% plot(vel,1-fittedint,'r')
% hold off
% xlabel('Radial velocity (km/s)')
% ylabel('Intensity')
% legend('Observed line','Fitted line','Location','SouthEast')
% axis([-inf inf -0.05 1.1])
% disp(' ')
% fprintf('\nResults:\nFWHM for line and spectrograph response fn = %.2f km/s\nV sin(i) = %.2f km/s\n\n',results1(1),results1(2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_vsini_directfit_inputs(vel,int,locs)

%Check velocity input
if ~isnumeric(vel) || ~isvector(vel), error('vel must be a numeric vector.'); end
if any(~isfinite(vel)), error('vel must contain only finite values.'); end

%Check intensity input
if ~isnumeric(int) || ~isvector(int), error('int must be a numeric vector.'); end
if any(~isfinite(int)), error('int must contain only finite values.'); end
if numel(int)~=numel(vel), error('vel and int must have the same number of elements.'); end

%Check line-centre input
if ~isnumeric(locs) || isempty(locs) || any(~isfinite(locs(:))), error('locs must contain at least one finite numeric line-centre value.'); end

end