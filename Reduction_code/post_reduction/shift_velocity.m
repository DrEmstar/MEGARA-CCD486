function [newwave]=shift_velocity(wave,vel)

%Validate inputs
if ~isnumeric(wave) || ~isvector(wave), error('wave must be a numeric vector.'); end
if any(~isfinite(wave(:))) || any(wave(:)<=0), error('wave must contain only finite positive values.'); end
if ~isnumeric(vel) || ~isscalar(vel) || ~isfinite(vel), error('vel must be a finite scalar numeric value.'); end

%Apply velocity shift
c=299792.458; %speed of light (km/s)
if vel<=-c, error('vel must be greater than -c.'); end
shift=log(1+vel/c); %compute shift in log space
newwave=exp(log(wave)-shift); %create new shifted axis in wavelength space

end
