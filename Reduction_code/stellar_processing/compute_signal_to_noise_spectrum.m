function signal_to_noise=compute_signal_to_noise_spectrum(stellar_order,exp_time)

%Validate inputs
validate_signal_to_noise_inputs(stellar_order,exp_time);

%Load photon conversion
photon_conversion=load_required_mat_file('photon_ADU_conversion.mat');
if ~isfield(photon_conversion,'mu_0_avg'), error('photon_ADU_conversion.mat must contain mu_0_avg.'); end
mu_0_avg=photon_conversion.mu_0_avg;

%dark noise
dark=find_dark_value(exp_time);

%Load readout/bias value
readout_noise_data=load_required_mat_file('readout_noise_value.mat');
if ~isfield(readout_noise_data,'bias'), error('readout_noise_value.mat must contain bias.'); end
bias=readout_noise_data.bias;

%Subtract bias contribution from dark estimate
dark=dark-bias;

%shotnoise- extracted order 100 from the stellar image
image_sum=mean(stellar_order.data,2,'omitnan');

%scaling backgrounds
detector_pixel_count=4096*4036;
bias=bias/detector_pixel_count;
sigma_bias=sqrt(max(bias,0));
dark=dark/detector_pixel_count;
sigma_dark=sqrt(max(dark,0));

%Calculate shot noise
shot_noise_argument=image_sum-bias-dark;
shot_noise_argument(shot_noise_argument<0)=0;
sigma_shot=sqrt(shot_noise_argument);

%backround total
B=sqrt(sigma_shot.^2+sigma_dark^2+sigma_bias^2)*mu_0_avg;

%light total
N=image_sum*mu_0_avg;
L=N-B;
L(L<0)=0;

%Calculate signal-to-noise spectrum
signal_to_noise_spectrum=L./sqrt(L+2*B);

%converting to a single value using the section of order 100 where there is no overlap
signal_to_noise_region_start=1700;
signal_to_noise_region_stop=2219;
if numel(signal_to_noise_spectrum)<signal_to_noise_region_start
    warning('Order is too short for the default signal-to-noise region. Using the full available order.')
    signal_to_noise_region=signal_to_noise_spectrum;
else
    signal_to_noise_region_stop=min(signal_to_noise_region_stop,numel(signal_to_noise_spectrum));
    signal_to_noise_region=signal_to_noise_spectrum(signal_to_noise_region_start:signal_to_noise_region_stop);
end

%Summarise signal-to-noise
signal_to_noise=median(signal_to_noise_region,'omitnan');
if ~isfinite(signal_to_noise), signal_to_noise=NaN; else, signal_to_noise=round(signal_to_noise,0); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_signal_to_noise_inputs(stellar_order,exp_time)

%Check stellar order
if ~isstruct(stellar_order), error('stellar_order must be a structure.'); end
if ~isfield(stellar_order,'data'), error('stellar_order must contain data.'); end
if ~isnumeric(stellar_order.data) || ndims(stellar_order.data)~=2, error('stellar_order.data must be a numeric 2D array.'); end
if isempty(stellar_order.data), error('stellar_order.data is empty.'); end

%Check exposure time
if ~isnumeric(exp_time) || ~isscalar(exp_time) || ~isfinite(exp_time) || exp_time<0, error('exp_time must be a finite scalar value >= 0.'); end

end

function loaded_file=load_required_mat_file(file_name)

%Find required MAT file in current folder or MATLAB path
if isfile(file_name)
    file_path=file_name;
else
    file_path=which(file_name);
end

%Stop if file cannot be found
if isempty(file_path)
    error(['Required MAT file not found on current folder or MATLAB path: ' file_name])
end

%Load required MAT file
loaded_file=load(file_path);

end