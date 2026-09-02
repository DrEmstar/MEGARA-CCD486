function dark=find_dark_value(exp_time)

%Validate input
if ~isnumeric(exp_time) || ~isscalar(exp_time) || ~isfinite(exp_time) || exp_time<0, error('exp_time must be a finite scalar value >= 0.'); end

%Find dark-exposure calibration file
dark_exposure_fit_file=which('dark_exposure_fit.mat');
if isempty(dark_exposure_fit_file), error('Could not find dark_exposure_fit.mat on the MATLAB path.'); end

%Load dark-exposure fit
loaded_dark_fit=load(dark_exposure_fit_file,'p','S','mu');
if ~isfield(loaded_dark_fit,'p') || ~isfield(loaded_dark_fit,'S') || ~isfield(loaded_dark_fit,'mu'), error('dark_exposure_fit.mat must contain p, S, and mu.'); end

%Calculate dark value
dark=polyval(loaded_dark_fit.p,exp_time,loaded_dark_fit.S,loaded_dark_fit.mu);

end