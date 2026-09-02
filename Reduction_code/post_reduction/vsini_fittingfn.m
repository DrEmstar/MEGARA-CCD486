function V=vsini_fittingfn(fit_parameters,xdata,lambda0,limb_darkening_coefficient)
%works out broadened line profile in lamda space and converts to velocity space
% x(1)=FWHM gaussian
% x(2)=V sin(i)

%Validate inputs
validate_vsini_fitting_inputs(fit_parameters,xdata,lambda0,limb_darkening_coefficient);

%Set constants
c=299792.458; %speed of light (km/s)

%Read fit parameters
gaussian_fwhm=fit_parameters(1);
projected_rotational_velocity=fit_parameters(2);

%Return invalid model for invalid trial values
if gaussian_fwhm<=0 || projected_rotational_velocity<=0
    V=nan(size(xdata));
    return
end

%Build fine wavelength grid
delta_lambda=-10:0.001:10; %make a fine wavelength grid on which to construct
delta_lambda_gaussian=-10:0.001:10; %our broadening curve

%Calculate rotational broadening kernel
rotational_velocity=(delta_lambda*c/lambda0)/projected_rotational_velocity;
valid_rotational_pixels=abs(rotational_velocity)<=1;
if ~any(valid_rotational_pixels)
    V=nan(size(xdata));
    return
end
delta_lambda_rotation=delta_lambda(valid_rotational_pixels);
rotational_velocity=rotational_velocity(valid_rotational_pixels);
rotational_kernel=(2*(1-limb_darkening_coefficient)*sqrt(1-rotational_velocity.^2)+0.5*pi*limb_darkening_coefficient*(1-rotational_velocity.^2))/(pi*lambda0*projected_rotational_velocity/c*(1-limb_darkening_coefficient/3));

%Calculate Gaussian kernel
gaussian_kernel=exp(-(2*sqrt(log(2))*delta_lambda_gaussian/gaussian_fwhm).^2);

%Convolve rotational and Gaussian kernels
broadened_kernel=conv(rotational_kernel,gaussian_kernel);

%Find convolution centre
[~,rotation_centre_index]=min(abs(delta_lambda_rotation));
[~,gaussian_centre_index]=min(abs(delta_lambda_gaussian));
convolution_centre_index=rotation_centre_index+gaussian_centre_index-1;

%make an x-axis for the convolution
number_left_pixels=convolution_centre_index-1;
number_right_pixels=numel(broadened_kernel)-convolution_centre_index;
new_wavelength_axis=(-number_left_pixels:number_right_pixels)*0.001;

%convert the x-axis to velocity space
new_velocity_axis=new_wavelength_axis/lambda0*c;

%project it on to the required axis
V=interp1(new_velocity_axis,broadened_kernel,xdata,'spline',0);

%normalise
maximum_model_value=max(V);
if ~isfinite(maximum_model_value) || maximum_model_value<=0
    V=nan(size(xdata));
else
    V=V/maximum_model_value; %normalise
end

% hold on
% plot(newax,V)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_vsini_fitting_inputs(fit_parameters,xdata,lambda0,limb_darkening_coefficient)

%Check fit parameters
if ~isnumeric(fit_parameters) || numel(fit_parameters)<2, error('fit_parameters must contain at least two numeric values.'); end
if any(~isfinite(fit_parameters(1:2))), error('fit_parameters must contain finite values.'); end

%Check velocity axis
if ~isnumeric(xdata) || ~isvector(xdata), error('xdata must be a numeric vector.'); end
if any(~isfinite(xdata)), error('xdata must contain only finite values.'); end

%Check wavelength
if ~isnumeric(lambda0) || ~isscalar(lambda0) || ~isfinite(lambda0) || lambda0<=0, error('lambda0 must be a finite positive scalar.'); end

%Check limb darkening
if ~isnumeric(limb_darkening_coefficient) || ~isscalar(limb_darkening_coefficient) || ~isfinite(limb_darkening_coefficient), error('limb_darkening_coefficient must be a finite scalar.'); end
if limb_darkening_coefficient<0 || limb_darkening_coefficient>1, error('limb_darkening_coefficient must be between 0 and 1.'); end

end