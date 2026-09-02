function [radial_velocity,rv_error,systemic_velocity_direct,systemic_velocity_gauss,vsini,vsini_error,vbroad,vbroad_error,vsini_fit]=radial_velocity_calculation(lsd_profile)

%RADIAL_VELOCITY_CALCULATION Measure RV and observed broadening from LSD profiles.
%VBROAD is the rotation-equivalent observed width, including instrumental
%broadening. VSINI has the calibrated instrumental variance removed, but
%macroturbulent broadening is not removed from either quantity.


%Standardise input profile
[velocity,intensity,jd]=standardise_lsd_profile(lsd_profile);
number_spectra=size(intensity,1);
velocity=velocity(:).';
velocity_step=median(abs(diff(velocity)),'omitnan');
if ~isfinite(velocity_step) || velocity_step<=0, velocity_step=1; end
[instrumental_sigma,instrumental_calibration_file]=load_instrumental_sigma();

%Initialise outputs
radial_velocity=nan(number_spectra,1);
rv_error=nan(number_spectra,1);
vsini=nan(number_spectra,1);
vsini_error=NaN;
vbroad=nan(number_spectra,1);
vbroad_error=NaN;
vsini_fit=struct();

%Build representative profile and common moment window
median_intensity=median(intensity,1,'omitnan');
median_absorption=make_absorption_profile(velocity,median_intensity,[]);
[line_window,representative_centre,representative_width]=find_common_line_window(velocity,median_absorption);
window_index=velocity>=line_window(1) & velocity<=line_window(2);
if sum(window_index)<8, error('Could not define a usable LSD line window for RV/vsini measurement.'); end
maximum_half_width=0.5*representative_width;
[integration_half_width,window_selection]=select_moment_half_width( ...
    velocity,median_intensity,representative_centre,maximum_half_width, ...
    instrumental_sigma);

%Measure the representative profile with the same moments used for every epoch
representative_moments=measure_profile_moments(velocity,median_intensity, ...
    representative_centre,integration_half_width,instrumental_sigma);
systemic_velocity_direct=representative_moments.first_moment;
systemic_velocity_gauss=measure_gaussian_centre(velocity,median_absorption, ...
    line_window,systemic_velocity_direct);
if ~isfinite(systemic_velocity_direct), error('Could not measure the median LSD first moment.'); end
if ~isfinite(systemic_velocity_gauss), systemic_velocity_gauss=systemic_velocity_direct; end
representative_vsini=representative_moments.equivalent_vsini;
representative_vbroad=representative_moments.equivalent_vbroad;

%Initialise diagnostic arrays
profile_centre=nan(number_spectra,1);
profile_centre_centroid=nan(number_spectra,1);
profile_depth=nan(number_spectra,1);
profile_fwhm=nan(number_spectra,1);
profile_percentile_width=nan(number_spectra,1);
fit_quality_flag=nan(number_spectra,1);
zeroth_moment=nan(number_spectra,1);
second_moment_observed=nan(number_spectra,1);
second_moment_stellar=nan(number_spectra,1);
third_moment=nan(number_spectra,1);
baseline_offset=nan(number_spectra,1);
baseline_slope=nan(number_spectra,1);
moment_window=nan(number_spectra,2);
vsini_error=nan(number_spectra,1);
vbroad_error=nan(number_spectra,1);

%Measure each profile independently with a moving moment window
for spectrum_index=1:number_spectra

    this_intensity=intensity(spectrum_index,:);
    this_moments=measure_profile_moments(velocity,this_intensity, ...
        systemic_velocity_direct,integration_half_width,instrumental_sigma);

    profile_centre(spectrum_index)=this_moments.first_moment;
    profile_centre_centroid(spectrum_index)=this_moments.first_moment;
    radial_velocity(spectrum_index)=this_moments.first_moment-systemic_velocity_direct;
    rv_error(spectrum_index)=this_moments.first_moment_error;
    vsini(spectrum_index)=this_moments.equivalent_vsini;
    vsini_error(spectrum_index)=this_moments.equivalent_vsini_error;
    vbroad(spectrum_index)=this_moments.equivalent_vbroad;
    vbroad_error(spectrum_index)=this_moments.equivalent_vbroad_error;
    fit_quality_flag(spectrum_index)=this_moments.quality_flag;
    profile_depth(spectrum_index)=this_moments.depth;
    profile_fwhm(spectrum_index)=this_moments.fwhm;
    profile_percentile_width(spectrum_index)=this_moments.percentile_width;
    zeroth_moment(spectrum_index)=this_moments.zeroth_moment;
    second_moment_observed(spectrum_index)=this_moments.second_moment_observed;
    second_moment_stellar(spectrum_index)=this_moments.second_moment_stellar;
    third_moment(spectrum_index)=this_moments.third_moment;
    baseline_offset(spectrum_index)=this_moments.baseline_offset;
    baseline_slope(spectrum_index)=this_moments.baseline_slope;
    moment_window(spectrum_index,:)=this_moments.window;

end

%Apply a sampling floor to formal first-moment errors
rv_error(~isfinite(rv_error) | rv_error<=0)=velocity_step;
rv_error=max(rv_error,0.05*velocity_step);

%Store diagnostics
vsini_fit.velocity=velocity;
vsini_fit.jd=jd;
vsini_fit.line_window=line_window;
vsini_fit.systemic_velocity_direct=systemic_velocity_direct;
vsini_fit.systemic_velocity_gauss=systemic_velocity_gauss;
vsini_fit.representative_vbroad=representative_vbroad;
vsini_fit.representative_vsini=representative_vsini;
vsini_fit.representative_moments=representative_moments;
vsini_fit.integration_half_width=integration_half_width;
vsini_fit.window_selection=window_selection;
vsini_fit.instrumental_sigma=instrumental_sigma;
vsini_fit.instrumental_fwhm=2.354820045*instrumental_sigma;
vsini_fit.instrumental_calibration_file=instrumental_calibration_file;
vsini_fit.profile_centre=profile_centre;
vsini_fit.profile_centre_centroid=profile_centre_centroid;
vsini_fit.profile_depth=profile_depth;
vsini_fit.profile_fwhm=profile_fwhm;
vsini_fit.profile_percentile_width=profile_percentile_width;
vsini_fit.vbroad=vbroad;
vsini_fit.vbroad_error=vbroad_error;
vsini_fit.vsini=vsini;
vsini_fit.vsini_error=vsini_error;
vsini_fit.fit_quality_flag=fit_quality_flag;
vsini_fit.zeroth_moment=zeroth_moment;
vsini_fit.second_moment_observed=second_moment_observed;
vsini_fit.second_moment_stellar=second_moment_stellar;
vsini_fit.third_moment=third_moment;
vsini_fit.baseline_offset=baseline_offset;
vsini_fit.baseline_slope=baseline_slope;
vsini_fit.moment_window=moment_window;
vsini_fit.method='iterative_first_and_second_moments';
vsini_fit.notes=['RV is the first moment relative to the median profile. ' ...
    'vbroad is the rotation-equivalent observed width and includes ' ...
    'instrumental broadening. vsini has the calibrated instrumental ' ...
    'variance removed. Macroturbulence is not removed from either value, ' ...
    'so vsini is not a pure projected rotational velocity.'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [velocity,intensity,jd]=standardise_lsd_profile(lsd_profile)

%Check fields
if ~isstruct(lsd_profile), error('lsd_profile must be a structure.'); end
if ~isfield(lsd_profile,'velocity') || ~isfield(lsd_profile,'intensity'), error('lsd_profile must contain velocity and intensity fields.'); end

velocity=lsd_profile.velocity(:).';
intensity=lsd_profile.intensity;

%Orient intensity as spectra x velocity
if size(intensity,2)==numel(velocity)
    %Already correct
elseif size(intensity,1)==numel(velocity)
    intensity=intensity.';
else
    error(['LSD intensity size is incompatible with velocity axis. Velocity has ' num2str(numel(velocity)) ' points; intensity is ' num2str(size(intensity,1)) ' x ' num2str(size(intensity,2)) '.'])
end

%Sort velocity axis
[velocity,sort_index]=sort(velocity);
intensity=intensity(:,sort_index);

%Read or create times
if isfield(lsd_profile,'jd')
    jd=lsd_profile.jd(:);
else
    jd=(1:size(intensity,1)).';
end

end

function moments=measure_profile_moments(velocity,intensity,centre_start,half_width,instrumental_sigma)

%Measure iterative moments after fitting a local linear wing baseline
velocity=velocity(:).';
intensity=intensity(:).';
finite_index=isfinite(velocity) & isfinite(intensity);
centre=centre_start;
maximum_iterations=8;
convergence_limit=0.01*median(abs(diff(velocity)),'omitnan');

for iteration=1:maximum_iterations
    integration_index=finite_index & abs(velocity-centre)<=half_width;
    wing_index=finite_index & abs(velocity-centre)>=1.10*half_width & ...
        abs(velocity-centre)<=1.60*half_width;
    if sum(wing_index)<8
        wing_index=finite_index & ~integration_index;
    end
    if sum(integration_index)<8 || sum(wing_index)<3
        moments=failed_moments(centre,half_width);
        return
    end

    wing_velocity=velocity(wing_index)-centre;
    baseline_parameters=[ones(sum(wing_index),1) wing_velocity(:)]\ ...
        intensity(wing_index).';
    baseline=baseline_parameters(1)+baseline_parameters(2).*(velocity-centre);
    absorption=baseline-intensity;
    delta_velocity=gradient(velocity);
    moment_weight=absorption(integration_index).*delta_velocity(integration_index);
    zeroth_moment=sum(moment_weight,'omitnan');
    if ~isfinite(zeroth_moment) || zeroth_moment<=0
        moments=failed_moments(centre,half_width);
        return
    end
    local_velocity=velocity(integration_index);
    new_centre=sum(local_velocity.*moment_weight,'omitnan')/zeroth_moment;
    if ~isfinite(new_centre) || abs(new_centre-centre)>half_width
        moments=failed_moments(centre,half_width);
        return
    end
    if abs(new_centre-centre)<=convergence_limit
        centre=new_centre;
        break
    end
    centre=new_centre;
end

%Recalculate the final moments in the converged moving window
integration_index=finite_index & abs(velocity-centre)<=half_width;
wing_index=finite_index & abs(velocity-centre)>=1.10*half_width & ...
    abs(velocity-centre)<=1.60*half_width;
if sum(wing_index)<8, wing_index=finite_index & ~integration_index; end
wing_velocity=velocity(wing_index)-centre;
baseline_parameters=[ones(sum(wing_index),1) wing_velocity(:)]\ ...
    intensity(wing_index).';
baseline=baseline_parameters(1)+baseline_parameters(2).*(velocity-centre);
absorption=baseline-intensity;
delta_velocity=gradient(velocity);
local_velocity=velocity(integration_index);
local_absorption=absorption(integration_index);
local_delta_velocity=delta_velocity(integration_index);
moment_weight=local_absorption.*local_delta_velocity;
zeroth_moment=sum(moment_weight,'omitnan');
first_moment=sum(local_velocity.*moment_weight,'omitnan')/zeroth_moment;
second_moment_observed=sum((local_velocity-first_moment).^2.* ...
    moment_weight,'omitnan')/zeroth_moment;
third_moment=sum((local_velocity-first_moment).^3.*moment_weight, ...
    'omitnan')/zeroth_moment;
second_moment_stellar=second_moment_observed-instrumental_sigma^2;

%Convert observed and instrument-corrected variances to equivalent rotation
%parameters for a linear limb-darkening coefficient epsilon=0.6.
rotation_variance_coefficient=0.225;
if isfinite(second_moment_observed) && second_moment_observed>0
    equivalent_vbroad=sqrt(second_moment_observed/rotation_variance_coefficient);
else
    equivalent_vbroad=0;
end
if isfinite(second_moment_stellar) && second_moment_stellar>0
    equivalent_vsini=sqrt(second_moment_stellar/rotation_variance_coefficient);
else
    equivalent_vsini=0;
end

%Propagate wing noise through the first and second moments
wing_residual=intensity(wing_index)-baseline(wing_index);
noise_sigma=robust_sigma(wing_residual);
if ~isfinite(noise_sigma) || noise_sigma<=0, noise_sigma=std(wing_residual,'omitnan'); end
first_variance=noise_sigma^2*sum(((local_velocity-first_moment).* ...
    local_delta_velocity).^2,'omitnan')/zeroth_moment^2;
second_variance=noise_sigma^2*sum((((local_velocity-first_moment).^2- ...
    second_moment_observed).*local_delta_velocity).^2,'omitnan')/ ...
    zeroth_moment^2;
first_moment_error=sqrt(max(0,first_variance));
if equivalent_vbroad>0 && isfinite(second_variance)
    equivalent_vbroad_error=sqrt(second_variance)/ ...
        (2*rotation_variance_coefficient*equivalent_vbroad);
else
    equivalent_vbroad_error=sqrt(max(0,second_variance)/ ...
        rotation_variance_coefficient);
end
if equivalent_vsini>0 && isfinite(second_variance)
    equivalent_vsini_error=sqrt(second_variance)/ ...
        (2*rotation_variance_coefficient*equivalent_vsini);
else
    equivalent_vsini_error=sqrt(max(0,second_variance)/ ...
        rotation_variance_coefficient);
end

%Profile and quality diagnostics
depth=max(local_absorption,[],'omitnan');
line_window=[centre-half_width centre+half_width];
fwhm=measure_fwhm(velocity,absorption,line_window);
percentile_width=measure_percentile_width(velocity,absorption,line_window,0.05,0.95);
edge_count=max(2,round(0.10*numel(local_absorption)));
edge_absorption=median([local_absorption(1:edge_count) ...
    local_absorption(end-edge_count+1:end)],'omitnan');
quality_flag=0;
if ~isfinite(depth) || depth<=3*noise_sigma, quality_flag=max(quality_flag,1); end
if isfinite(depth) && edge_absorption>0.10*depth, quality_flag=max(quality_flag,2); end
if ~isfinite(second_moment_observed) || second_moment_observed<=0, quality_flag=3; end

moments=struct();
moments.zeroth_moment=zeroth_moment;
moments.first_moment=first_moment;
moments.first_moment_error=first_moment_error;
moments.second_moment_observed=second_moment_observed;
moments.second_moment_stellar=max(0,second_moment_stellar);
moments.third_moment=third_moment;
moments.equivalent_vbroad=equivalent_vbroad;
moments.equivalent_vbroad_error=equivalent_vbroad_error;
moments.equivalent_vsini=equivalent_vsini;
moments.equivalent_vsini_error=equivalent_vsini_error;
moments.baseline_offset=baseline_parameters(1);
moments.baseline_slope=baseline_parameters(2);
moments.depth=depth;
moments.fwhm=fwhm;
moments.percentile_width=percentile_width;
moments.window=[centre-half_width centre+half_width];
moments.quality_flag=quality_flag;
moments.noise_sigma=noise_sigma;
moments.edge_fraction=edge_absorption/depth;

end

function moments=failed_moments(centre,half_width)

%Return a consistent failed-moment result
moments=struct('zeroth_moment',NaN,'first_moment',NaN, ...
    'first_moment_error',NaN,'second_moment_observed',NaN, ...
    'second_moment_stellar',NaN,'third_moment',NaN, ...
    'equivalent_vbroad',NaN,'equivalent_vbroad_error',NaN, ...
    'equivalent_vsini',NaN,'equivalent_vsini_error',NaN, ...
    'baseline_offset',NaN,'baseline_slope',NaN,'depth',NaN, ...
    'fwhm',NaN,'percentile_width',NaN, ...
    'window',[centre-half_width centre+half_width], ...
    'quality_flag',9,'noise_sigma',NaN,'edge_fraction',NaN);

end

function [selected_half_width,selection]=select_moment_half_width( ...
    velocity,intensity,centre_start,maximum_half_width,instrumental_sigma)

%Select the smallest profile-contained window and add one grid-step margin
velocity_step=median(abs(diff(velocity)),'omitnan');
width_step=max(2,2*velocity_step);
minimum_half_width=max(10,8*velocity_step);
if maximum_half_width<minimum_half_width
    candidate_half_widths=maximum_half_width;
else
    candidate_half_widths=minimum_half_width:width_step:maximum_half_width;
    if candidate_half_widths(end)<maximum_half_width
        candidate_half_widths(end+1)=maximum_half_width;
    end
end

edge_fraction=nan(size(candidate_half_widths));
equivalent_width=nan(size(candidate_half_widths));
first_moment=nan(size(candidate_half_widths));
for width_index=1:numel(candidate_half_widths)
    this_moments=measure_profile_moments(velocity,intensity,centre_start, ...
        candidate_half_widths(width_index),instrumental_sigma);
    edge_fraction(width_index)=this_moments.edge_fraction;
    equivalent_width(width_index)=this_moments.zeroth_moment;
    first_moment(width_index)=this_moments.first_moment;
end

contained_index=find(isfinite(edge_fraction) & edge_fraction<=0.05 & ...
    isfinite(equivalent_width) & equivalent_width>0,1,'first');
if isempty(contained_index)
    selected_index=numel(candidate_half_widths);
    selection_status='No 5-percent edge crossing; used full detected line window';
else
    selected_index=min(contained_index+1,numel(candidate_half_widths));
    selection_status='First 5-percent edge crossing plus one width-grid safety step';
end
selected_half_width=candidate_half_widths(selected_index);

selection=struct();
selection.candidate_half_widths=candidate_half_widths;
selection.edge_fraction=edge_fraction;
selection.equivalent_width=equivalent_width;
selection.first_moment=first_moment;
selection.selected_index=selected_index;
selection.selected_half_width=selected_half_width;
selection.edge_limit=0.05;
selection.width_step=width_step;
selection.status=selection_status;

end

function absorption=make_absorption_profile(velocity,intensity,line_window)

%Estimate continuum from wings and convert an absorption line to a positive profile
velocity=velocity(:).';
intensity=intensity(:).';
finite_index=isfinite(velocity) & isfinite(intensity);

if nargin<3 || isempty(line_window)
    number_points=numel(velocity);
    wing_count=max(3,round(0.15*number_points));
    wing_index=false(size(velocity));
    wing_index(1:wing_count)=true;
    wing_index(end-wing_count+1:end)=true;
else
    line_width=line_window(2)-line_window(1);
    wing_index=finite_index & (velocity<line_window(1)-0.15*line_width | velocity>line_window(2)+0.15*line_width);
    if sum(wing_index)<6
        wing_index=finite_index & (velocity<line_window(1) | velocity>line_window(2));
    end
end

if sum(wing_index)>=3
    continuum=median(intensity(wing_index),'omitnan');
else
    continuum=max(intensity(finite_index),[],'omitnan');
end

absorption=continuum-intensity;
absorption(~finite_index)=NaN;

%Suppress tiny negative wing excursions without clipping real profile structure
negative_floor=-0.20*max(absorption,[],'omitnan');
if isfinite(negative_floor)
    absorption(absorption<negative_floor)=negative_floor;
end

end

function [line_window,representative_centre,representative_width]=find_common_line_window(velocity,absorption)

%Find a common line window from the representative absorption profile
finite_index=isfinite(velocity) & isfinite(absorption);
if sum(finite_index)<8, error('Cannot determine line window from fewer than eight finite points.'); end

velocity=velocity(:).';
absorption=absorption(:).';
smoothed_absorption=smooth_vector(absorption,max(3,round(numel(absorption)/100)));
[peak_absorption,peak_index]=max(smoothed_absorption,[],'omitnan');

if ~isfinite(peak_absorption) || peak_absorption<=0
    representative_centre=median(velocity(finite_index),'omitnan');
    representative_width=0.5*(max(velocity(finite_index))-min(velocity(finite_index)));
    line_window=[representative_centre-representative_width representative_centre+representative_width];
    return
end

representative_centre=velocity(peak_index);
threshold=max(0.05*peak_absorption,median(smoothed_absorption(finite_index),'omitnan')+2*robust_sigma(smoothed_absorption(finite_index)));
above_threshold=smoothed_absorption>threshold & finite_index;

%Keep contiguous region around the peak
left_index=peak_index;
while left_index>1 && above_threshold(left_index-1)
    left_index=left_index-1;
end
right_index=peak_index;
while right_index<numel(velocity) && above_threshold(right_index+1)
    right_index=right_index+1;
end

%Expand the window to include the full line wings
velocity_step=median(abs(diff(velocity)),'omitnan');
if ~isfinite(velocity_step) || velocity_step<=0, velocity_step=1; end
expand_pixels=max(3,round(10/velocity_step));
left_index=max(1,left_index-expand_pixels);
right_index=min(numel(velocity),right_index+expand_pixels);

line_window=[velocity(left_index) velocity(right_index)];
representative_width=line_window(2)-line_window(1);

%Fallback if window is too small
if representative_width<8*velocity_step
    profile_width=measure_percentile_width(velocity,absorption,[min(velocity) max(velocity)],0.05,0.95);
    if isfinite(profile_width) && profile_width>8*velocity_step
        representative_width=profile_width;
        line_window=[representative_centre-0.65*profile_width representative_centre+0.65*profile_width];
    else
        representative_width=max(40,20*velocity_step);
        line_window=[representative_centre-representative_width/2 representative_centre+representative_width/2];
    end
end

line_window(1)=max(line_window(1),min(velocity));
line_window(2)=min(line_window(2),max(velocity));

end

function centre=measure_absorption_centroid(velocity,absorption,line_window,default_centre)

%Measure absorption-weighted centroid in the line window
window_index=velocity>=line_window(1) & velocity<=line_window(2) & isfinite(absorption);
if sum(window_index)<3, centre=default_centre; return; end
local_velocity=velocity(window_index);
local_absorption=absorption(window_index);
local_absorption=local_absorption-min(0,min(local_absorption,[],'omitnan'));
local_absorption(local_absorption<0)=0;
if sum(local_absorption,'omitnan')<=0
    centre=default_centre;
else
    centre=sum(local_velocity.*local_absorption,'omitnan')./sum(local_absorption,'omitnan');
end

end

function centre=measure_gaussian_centre(velocity,absorption,line_window,centre_start)

%Fit a simple Gaussian absorption profile to estimate centre only
window_index=velocity>=line_window(1) & velocity<=line_window(2) & isfinite(absorption);
local_velocity=velocity(window_index);
local_absorption=absorption(window_index);

if numel(local_velocity)<6 || all(~isfinite(local_absorption))
    centre=centre_start;
    return
end

local_absorption=local_absorption(:).';
local_velocity=local_velocity(:).';
peak_depth=max(local_absorption,[],'omitnan');
if ~isfinite(peak_depth) || peak_depth<=0
    centre=centre_start;
    return
end

sigma_start=max(2*median(abs(diff(local_velocity)),'omitnan'),0.25*(line_window(2)-line_window(1)));
start_parameters=[centre_start log(sigma_start) peak_depth median(local_absorption,'omitnan')];
objective=@(parameters) sum((local_absorption-gaussian_absorption_model(local_velocity,parameters)).^2,'omitnan');
options=optimset('Display','off','MaxIter',400,'MaxFunEvals',1000);
fit_parameters=fminsearch(objective,start_parameters,options);
centre=fit_parameters(1);

if ~isfinite(centre) || centre<line_window(1) || centre>line_window(2)
    centre=centre_start;
end

end

function model=gaussian_absorption_model(velocity,parameters)

%Gaussian absorption model
centre=parameters(1);
sigma=exp(parameters(2));
depth=parameters(3);
offset=parameters(4);
model=offset+depth.*exp(-0.5*((velocity-centre)./sigma).^2);

end

function fit_result=fit_rotational_profile(velocity,absorption,line_window,centre_start,vsini_start,velocity_step,vsini_bounds,instrumental_sigma)

%Fit a rotational profile to one absorption profile using multiple starts
velocity=velocity(:).';
absorption=absorption(:).';
window_index=velocity>=line_window(1) & velocity<=line_window(2) & isfinite(absorption);
local_velocity=velocity(window_index);
local_absorption=absorption(window_index);

if numel(local_velocity)<8 || all(~isfinite(local_absorption))
    fit_result=make_failed_fit_result(velocity,centre_start,vsini_start);
    return
end

%Clean local profile
local_absorption=local_absorption(:).';
local_velocity=local_velocity(:).';
local_depth=max(local_absorption,[],'omitnan')-median(local_absorption,'omitnan');
if ~isfinite(local_depth) || local_depth<=0, local_depth=max(local_absorption,[],'omitnan'); end
if ~isfinite(local_depth) || local_depth<=0
    fit_result=make_failed_fit_result(velocity,centre_start,vsini_start);
    return
end

%Define broad data-derived bounds
line_width=line_window(2)-line_window(1);
if nargin<7 || isempty(vsini_bounds)
    lower_vsini=0;
    upper_vsini=max(lower_vsini+velocity_step,1.2*line_width);
else
    lower_vsini=vsini_bounds(1);
    upper_vsini=vsini_bounds(2);
end

lower_centre=line_window(1)+2*velocity_step;
upper_centre=line_window(2)-2*velocity_step;
centre_start=min(max(centre_start,lower_centre),upper_centre);
vsini_start=min(max(vsini_start,lower_vsini),upper_vsini);

%Profile weights: include the full line but mildly downweight the far continuum
normalised_absorption=local_absorption./max(local_absorption,[],'omitnan');
profile_weights=0.35+0.65*max(0,normalised_absorption);
profile_weights(~isfinite(profile_weights))=1;

%Multiple starting points protect against local minima
centre_offsets=[0 -2*velocity_step 2*velocity_step -5*velocity_step 5*velocity_step];
vsini_factors=[0.65 0.85 1.0 1.25 1.60];

best_objective=Inf;
best_parameters=[centre_start vsini_start];
best_linear_parameters=[local_depth 0 0];
quality_flag=0;
optimisation_options=optimset('Display','off','MaxIter',500,'MaxFunEvals',1500);

for centre_index=1:numel(centre_offsets)
    for vsini_index=1:numel(vsini_factors)
        trial_centre=min(max(centre_start+centre_offsets(centre_index),lower_centre),upper_centre);
        trial_vsini=min(max(vsini_start*vsini_factors(vsini_index),lower_vsini),upper_vsini);
        trial_parameters=[trial_centre trial_vsini];
        objective=@(parameters) rotational_profile_objective(parameters,local_velocity,local_absorption,profile_weights,lower_centre,upper_centre,lower_vsini,upper_vsini,instrumental_sigma);
        [fit_parameters,fit_objective]=fminsearch(objective,trial_parameters,optimisation_options);
        [fit_objective,linear_parameters]=rotational_profile_objective(fit_parameters,local_velocity,local_absorption,profile_weights,lower_centre,upper_centre,lower_vsini,upper_vsini,instrumental_sigma);
        if fit_objective<best_objective
            best_objective=fit_objective;
            best_parameters=fit_parameters;
            best_linear_parameters=linear_parameters;
        end
    end
end

%Apply bounds exactly after optimisation
centre=min(max(best_parameters(1),lower_centre),upper_centre);
fit_vsini=min(max(abs(best_parameters(2)),lower_vsini),upper_vsini);
[kernel,full_model_absorption]=make_full_rotational_model(velocity,centre,fit_vsini,best_linear_parameters,instrumental_sigma);
[~,local_model_absorption]=make_full_rotational_model(local_velocity,centre,fit_vsini,best_linear_parameters,instrumental_sigma);
residual=local_absorption-local_model_absorption;
chi_square=sqrt(mean(residual.^2,'omitnan'));
centre_error=estimate_centre_error(local_velocity,local_absorption,local_model_absorption,velocity_step);

%Quality flag: 0 good, 1 at/beside bound, 2 weak/shallow line, 3 poor residual
if fit_vsini<=lower_vsini+0.5*velocity_step || fit_vsini>=upper_vsini-0.5*velocity_step
    quality_flag=max(quality_flag,1);
end
if best_linear_parameters(1)<=0 || local_depth<=0
    quality_flag=max(quality_flag,2);
end
if isfinite(chi_square) && isfinite(local_depth) && chi_square>0.35*local_depth
    quality_flag=max(quality_flag,3);
end

fit_result=struct();
fit_result.centre=centre;
fit_result.vsini=fit_vsini;
fit_result.depth=best_linear_parameters(1);
fit_result.offset=best_linear_parameters(2);
fit_result.slope=best_linear_parameters(3);
fit_result.kernel=kernel;
fit_result.model_absorption=full_model_absorption;
fit_result.chi_square=chi_square;
fit_result.centre_error=centre_error;
fit_result.quality_flag=quality_flag;

end

function [objective_value,linear_parameters]=rotational_profile_objective(parameters,velocity,absorption,weights,lower_centre,upper_centre,lower_vsini,upper_vsini,instrumental_sigma)

%Objective for rotational profile fit with linear depth/offset/slope solved analytically
centre=parameters(1);
vsini=abs(parameters(2));
penalty=0;

if centre<lower_centre, penalty=penalty+1e6*(lower_centre-centre)^2; centre=lower_centre; end
if centre>upper_centre, penalty=penalty+1e6*(centre-upper_centre)^2; centre=upper_centre; end
if vsini<lower_vsini, penalty=penalty+1e6*(lower_vsini-vsini)^2; vsini=lower_vsini; end
if vsini>upper_vsini, penalty=penalty+1e6*(vsini-upper_vsini)^2; vsini=upper_vsini; end

kernel=instrument_convolved_rotational_kernel(velocity-centre,vsini,instrumental_sigma);
velocity_scaled=(velocity-centre)./max(1,vsini);
design_matrix=[kernel(:) ones(numel(velocity),1) velocity_scaled(:)];
absorption_vector=absorption(:);
weight_vector=weights(:);
valid_index=isfinite(absorption_vector) & isfinite(weight_vector);
weighted_design=design_matrix(valid_index,:).*weight_vector(valid_index);
weighted_absorption=absorption_vector(valid_index).*weight_vector(valid_index);

if size(weighted_design,1)<3
    linear_parameters=[NaN NaN NaN];
    objective_value=Inf;
    return
end

linear_parameters=weighted_design\weighted_absorption;

%Prevent negative absorption depth from being selected
if linear_parameters(1)<0
    penalty=penalty+1e6*linear_parameters(1)^2;
end

model=design_matrix*linear_parameters;
residual=absorption_vector-model;
robust_scale=robust_sigma(residual(valid_index));
if ~isfinite(robust_scale) || robust_scale<=0, robust_scale=std(residual(valid_index),'omitnan'); end
if ~isfinite(robust_scale) || robust_scale<=0, robust_scale=1; end
huber_residual=huber_weighted_residual(residual(valid_index)./robust_scale,1.5).*robust_scale;
objective_value=sum((huber_residual.*weight_vector(valid_index)).^2,'omitnan')+penalty;
linear_parameters=linear_parameters(:).';

end

function [kernel,model_absorption]=make_full_rotational_model(velocity,centre,vsini,linear_parameters,instrumental_sigma)

%Evaluate rotational model on a velocity grid
kernel=instrument_convolved_rotational_kernel(velocity-centre,vsini,instrumental_sigma);
velocity_scaled=(velocity-centre)./max(1,vsini);
model_absorption=linear_parameters(1).*kernel+linear_parameters(2)+linear_parameters(3).*velocity_scaled;

end

function kernel=instrument_convolved_rotational_kernel(delta_velocity,vsini,instrumental_sigma)

%Convolve the stellar broadening kernel with the measured instrument profile
kernel=rotational_kernel_sampled(delta_velocity,vsini);
velocity_step=median(abs(diff(delta_velocity)),'omitnan');
if ~isfinite(instrumental_sigma) || instrumental_sigma<=0 || ...
        ~isfinite(velocity_step) || velocity_step<=0
    return
end
half_width=max(3,ceil(5*instrumental_sigma/velocity_step));
instrument_velocity=(-half_width:half_width).*velocity_step;
instrument_kernel=exp(-0.5*(instrument_velocity./instrumental_sigma).^2);
instrument_kernel=instrument_kernel./sum(instrument_kernel,'omitnan');
kernel=conv(kernel,instrument_kernel,'same');
normalisation=max(kernel,[],'omitnan');
if isfinite(normalisation) && normalisation>0, kernel=kernel./normalisation; end

end

function kernel=rotational_kernel_sampled(delta_velocity,vsini)

%Rotational broadening kernel sampled on the velocity grid and normalised to unit peak
limb_darkening=0.6;
velocity_step=median(abs(diff(delta_velocity)),'omitnan');
if vsini<0.25*velocity_step
    kernel=zeros(size(delta_velocity));
    [~,zero_index]=min(abs(delta_velocity));
    kernel(zero_index)=1;
    return
end
x=delta_velocity./max(vsini,eps);
kernel=zeros(size(delta_velocity));
inside=abs(x)<1;
mu=sqrt(max(0,1-x(inside).^2));
kernel(inside)=2*(1-limb_darkening).*mu+0.5*pi*limb_darkening.*mu.^2;
normalisation=max(kernel,[],'omitnan');
if isfinite(normalisation) && normalisation>0
    kernel=kernel./normalisation;
end

end

function [instrumental_sigma,calibration_file]=load_instrumental_sigma()

%Load the compact release calibration measured from extracted ThAr lines.
calibration_file=fullfile(fileparts(mfilename('fullpath')), ...
    'instrumental_broadening_calibration.mat');
if ~isfile(calibration_file)
    error(['Instrumental broadening calibration not found: ' calibration_file])
end
loaded_calibration=load(calibration_file,'instrumental_broadening');
if ~isfield(loaded_calibration,'instrumental_broadening') || ...
        ~isfield(loaded_calibration.instrumental_broadening,'sigma_velocity')
    error('instrumental_broadening_calibration.mat is missing sigma_velocity.')
end
instrumental_sigma=loaded_calibration.instrumental_broadening.sigma_velocity;
if ~isnumeric(instrumental_sigma) || ~isscalar(instrumental_sigma) || ...
        ~isfinite(instrumental_sigma) || instrumental_sigma<=0
    error('The instrumental sigma must be a finite scalar value > 0.')
end

end

function fit_result=make_failed_fit_result(velocity,centre_start,vsini_start)

%Return a failed-fit structure with consistent fields
fit_result=struct();
fit_result.centre=centre_start;
fit_result.vsini=vsini_start;
fit_result.depth=NaN;
fit_result.offset=NaN;
fit_result.slope=NaN;
fit_result.kernel=nan(size(velocity));
fit_result.model_absorption=nan(size(velocity));
fit_result.chi_square=NaN;
fit_result.centre_error=NaN;
fit_result.quality_flag=9;

end

function vsini_estimate=estimate_vsini_from_width(velocity,absorption,line_window,representative_width)

%Estimate initial vsini from robust full-profile widths
fwhm=measure_fwhm(velocity,absorption,line_window);
percentile_width=measure_percentile_width(velocity,absorption,line_window,0.05,0.95);
width_candidates=[fwhm percentile_width/1.6 representative_width/1.5];
width_candidates=width_candidates(isfinite(width_candidates) & width_candidates>0);
if isempty(width_candidates)
    vsini_estimate=NaN;
else
    vsini_estimate=0.5*median(width_candidates,'omitnan');
end

end

function width=measure_fwhm(velocity,absorption,line_window)

%Measure full width at half maximum for positive absorption profile
window_index=velocity>=line_window(1) & velocity<=line_window(2) & isfinite(absorption);
local_velocity=velocity(window_index);
local_absorption=absorption(window_index);
if numel(local_velocity)<5, width=NaN; return; end
local_absorption=local_absorption-min(local_absorption,[],'omitnan');
peak=max(local_absorption,[],'omitnan');
if ~isfinite(peak) || peak<=0, width=NaN; return; end
half_level=0.5*peak;
above=local_absorption>=half_level;
if sum(above)<2, width=NaN; return; end
first_index=find(above,1,'first');
last_index=find(above,1,'last');
left_velocity=interpolate_crossing(local_velocity,local_absorption,first_index,half_level,'left');
right_velocity=interpolate_crossing(local_velocity,local_absorption,last_index,half_level,'right');
width=right_velocity-left_velocity;
if width<=0, width=NaN; end

end

function crossing_velocity=interpolate_crossing(velocity,absorption,index_value,level,side)

%Interpolate threshold crossing
if strcmp(side,'left')
    if index_value<=1, crossing_velocity=velocity(index_value); return; end
    x1=velocity(index_value-1); x2=velocity(index_value); y1=absorption(index_value-1); y2=absorption(index_value);
else
    if index_value>=numel(velocity), crossing_velocity=velocity(index_value); return; end
    x1=velocity(index_value); x2=velocity(index_value+1); y1=absorption(index_value); y2=absorption(index_value+1);
end
if y2==y1
    crossing_velocity=0.5*(x1+x2);
else
    crossing_velocity=x1+(level-y1).*(x2-x1)./(y2-y1);
end

end

function width=measure_percentile_width(velocity,absorption,line_window,lower_fraction,upper_fraction)

%Measure width enclosing a chosen fraction of the absorption area
window_index=velocity>=line_window(1) & velocity<=line_window(2) & isfinite(absorption);
local_velocity=velocity(window_index);
local_absorption=absorption(window_index);
if numel(local_velocity)<5, width=NaN; return; end
local_absorption=local_absorption-min(0,min(local_absorption,[],'omitnan'));
local_absorption(local_absorption<0)=0;
if sum(local_absorption,'omitnan')<=0, width=NaN; return; end
cumulative_absorption=cumsum(local_absorption,'omitnan');
cumulative_absorption=cumulative_absorption./cumulative_absorption(end);
[unique_cumulative,unique_index]=unique(cumulative_absorption,'stable');
unique_velocity=local_velocity(unique_index);
if numel(unique_cumulative)<2, width=NaN; return; end
lower_velocity=interp1(unique_cumulative,unique_velocity,lower_fraction,'linear','extrap');
upper_velocity=interp1(unique_cumulative,unique_velocity,upper_fraction,'linear','extrap');
width=upper_velocity-lower_velocity;
if width<=0, width=NaN; end

end

function centre_error=estimate_centre_error(velocity,absorption,model_absorption,velocity_step)

%Estimate practical centre uncertainty from residual scale and line gradient
residual=absorption-model_absorption;
residual_scale=robust_sigma(residual(isfinite(residual)));
line_gradient=gradient(model_absorption,velocity);
gradient_scale=max(abs(line_gradient),[],'omitnan');
if ~isfinite(residual_scale) || ~isfinite(gradient_scale) || gradient_scale<=0
    centre_error=velocity_step;
else
    centre_error=max(0.25*velocity_step,min(10*velocity_step,residual_scale./gradient_scale));
end

end

function output=smooth_vector(input,half_width)

%Smooth a vector with a moving mean while preserving size
if nargin<2 || half_width<1, output=input; return; end
window_size=max(3,2*round(half_width)+1);
output=movmean(input,window_size,'omitnan');

end

function sigma=robust_sigma(values)

%Robust standard-deviation estimate
values=values(:);
values=values(isfinite(values));
if isempty(values)
    sigma=NaN;
else
    median_value=median(values,'omitnan');
    sigma=1.4826*median(abs(values-median_value),'omitnan');
    if ~isfinite(sigma) || sigma==0
        sigma=std(values,'omitnan');
    end
end

end

function residual=huber_weighted_residual(scaled_residual,tuning_constant)

%Huber residual for robust least squares
residual=scaled_residual;
large=abs(scaled_residual)>tuning_constant;
residual(large)=sign(scaled_residual(large)).*sqrt(tuning_constant.*abs(scaled_residual(large)));

end
