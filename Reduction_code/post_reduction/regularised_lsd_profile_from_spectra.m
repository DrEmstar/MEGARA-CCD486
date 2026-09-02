function regularised_lsd_profile=regularised_lsd_profile_from_spectra(fullwave,fullint,lsd_line_mask,jd,regularised_lsd_options)

%Create regularised LSD profiles by solving a linear inverse problem
%
% Model:
%   absorption spectrum = line-mask design matrix * mean absorption profile
%
% The saved intensity convention is:
%   intensity = 1 - fitted_absorption_profile
%
% This gives a pseudo-line profile with continuum near 1 and absorption as a dip.

%Constants
speed_of_light=299792.458; % km/s

%Standardise inputs
[fullwave,fullint]=standardise_spectral_arrays(fullwave,fullint);
jd=jd(:);

if numel(jd)~=size(fullint,1)
    warning('Number of JD values does not match number of spectra. Replacing JD with NaN values.')
    jd=nan(size(fullint,1),1);
end

%Read options
velocity_limit=get_option_value(regularised_lsd_options,'velocity_limit',200);
velocity_step=get_option_value(regularised_lsd_options,'velocity_step',1.0);
regularisation_lambda=get_option_value(regularised_lsd_options,'regularisation_lambda',1e-2);
regularisation_order=get_option_value(regularised_lsd_options,'regularisation_order',2);
normalisation=get_option_value(regularised_lsd_options,'normalisation','unit_depth');
profile_baseline=get_option_value(regularised_lsd_options,'profile_baseline','zero');

%Create velocity grid
velocity=-velocity_limit:velocity_step:velocity_limit;
number_velocity_points=numel(velocity);

%Extract line mask
[line_wavelength,line_weight,reference_weight]=validate_lsd_line_mask(lsd_line_mask);

%Build sparse design matrix
[design_matrix,used_pixel_index,normalised_line_weight]=build_regularised_lsd_design_matrix(fullwave,line_wavelength,line_weight,velocity,speed_of_light,reference_weight);

if isempty(used_pixel_index)
    error('No spectral pixels overlap the LSD line mask within the requested velocity range.')
end

%Build regularisation matrix
regularisation_matrix=make_regularisation_matrix(number_velocity_points,regularisation_order);

%Build normal equations
left_hand_side=design_matrix.'*design_matrix+regularisation_lambda*(regularisation_matrix.'*regularisation_matrix)+eps*speye(number_velocity_points);

%Try to create a reusable solver
try
    regularised_solver=decomposition(left_hand_side,'chol');
    use_decomposition=true;
catch
    regularised_solver=[];
    use_decomposition=false;
end

%Solve one LSD profile per spectrum
number_spectra=size(fullint,1);
absorption_profile=nan(number_spectra,number_velocity_points);

for spectrum_index=1:number_spectra

    %Use absorption rather than intensity
    observed_absorption=1-fullint(spectrum_index,used_pixel_index).';
    finite_pixels=isfinite(observed_absorption);

    if sum(finite_pixels)<number_velocity_points
        warning(['Spectrum ' num2str(spectrum_index) ' has too few finite pixels for regularised LSD.'])
        continue
    end

    if all(finite_pixels)
        right_hand_side=design_matrix.'*observed_absorption;
        if use_decomposition
            fitted_profile=regularised_solver\right_hand_side;
        else
            fitted_profile=left_hand_side\right_hand_side;
        end
    else
        local_design_matrix=design_matrix(finite_pixels,:);
        local_regularised_matrix=local_design_matrix.'*local_design_matrix+regularisation_lambda*(regularisation_matrix.'*regularisation_matrix)+eps*speye(number_velocity_points);
        right_hand_side=local_design_matrix.'*observed_absorption(finite_pixels);
        fitted_profile=local_regularised_matrix\right_hand_side;
    end

    absorption_profile(spectrum_index,:)=fitted_profile(:).';

end

%Apply baseline convention
[absorption_profile,baseline_values]=apply_lsd_baseline(velocity,absorption_profile,profile_baseline);

%Apply normalisation convention
[absorption_profile,normalisation_scale]=normalise_lsd_absorption_profiles(velocity,absorption_profile,normalisation);

%Convert absorption profile to pseudo-intensity profile
intensity_profile=1-absorption_profile;

%Create output structure
regularised_lsd_profile=struct();
regularised_lsd_profile.velocity=velocity;
regularised_lsd_profile.intensity=intensity_profile;
regularised_lsd_profile.absorption=absorption_profile;
regularised_lsd_profile.jd=jd;

%Store method metadata
regularised_lsd_profile.profile_type='LSD';
regularised_lsd_profile.method='regularised_inverse_problem';
regularised_lsd_profile.model='absorption_spectrum = design_matrix * mean_absorption_profile';
regularised_lsd_profile.profile_convention='intensity = 1 - fitted_absorption_profile';
regularised_lsd_profile.velocity_convention='c * log(lambda_observed/lambda_rest)';
regularised_lsd_profile.normalisation=normalisation;
regularised_lsd_profile.profile_baseline=profile_baseline;
regularised_lsd_profile.normalisation_scale=normalisation_scale;
regularised_lsd_profile.baseline_values=baseline_values;

%Store regularisation metadata
regularised_lsd_profile.regularisation_lambda=regularisation_lambda;
regularised_lsd_profile.regularisation_order=regularisation_order;
regularised_lsd_profile.velocity_limit=velocity_limit;
regularised_lsd_profile.velocity_step=velocity_step;

%Store mask metadata
regularised_lsd_profile.line_mask_wavelength=line_wavelength;
regularised_lsd_profile.line_mask_weight=line_weight;
regularised_lsd_profile.line_mask_weight_normalised=normalised_line_weight;
regularised_lsd_profile.line_mask_reference_weight=reference_weight;
regularised_lsd_profile.number_mask_lines=numel(line_wavelength);
regularised_lsd_profile.number_spectral_pixels=numel(fullwave);
regularised_lsd_profile.number_used_spectral_pixels=numel(used_pixel_index);
regularised_lsd_profile.used_pixel_index=used_pixel_index;

if isfield(lsd_line_mask,'method'), regularised_lsd_profile.line_mask_method=lsd_line_mask.method; end
if isfield(lsd_line_mask,'weak_limit'), regularised_lsd_profile.line_mask_weak_limit=lsd_line_mask.weak_limit; end
if isfield(lsd_line_mask,'region_limits'), regularised_lsd_profile.line_mask_region_limits=lsd_line_mask.region_limits; end
if isfield(lsd_line_mask,'region_labels'), regularised_lsd_profile.line_mask_region_labels=lsd_line_mask.region_labels; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fullwave,fullint]=standardise_spectral_arrays(fullwave,fullint)

%Force wavelength vector into row orientation
fullwave=fullwave(:).';

%Check intensity matrix orientation
if size(fullint,2)==numel(fullwave)
    %Already spectra x wavelength
elseif size(fullint,1)==numel(fullwave)
    fullint=fullint.';
else
    error(['fullwave/fullint size mismatch. fullwave has ' num2str(numel(fullwave)) ' pixels; fullint is ' num2str(size(fullint,1)) ' x ' num2str(size(fullint,2)) '.'])
end

%Sort by wavelength if needed
[fullwave,sort_index]=sort(fullwave);
fullint=fullint(:,sort_index);

end

function value=get_option_value(options,field_name,default_value)

%Read option field or default
if isstruct(options) && isfield(options,field_name) && ~isempty(options.(field_name))
    value=options.(field_name);
else
    value=default_value;
end

end

function [line_wavelength,line_weight,reference_weight]=validate_lsd_line_mask(lsd_line_mask)

%Extract and validate LSD line mask
if ~isstruct(lsd_line_mask)
    error('lsd_line_mask must be a structure.')
end

if ~isfield(lsd_line_mask,'wavelength')
    error('lsd_line_mask is missing required field wavelength.')
end

if isfield(lsd_line_mask,'weight')
    line_weight=lsd_line_mask.weight(:);
elseif isfield(lsd_line_mask,'width')
    line_weight=lsd_line_mask.width(:);
else
    error('lsd_line_mask is missing required field weight or width.')
end

line_wavelength=lsd_line_mask.wavelength(:);

if numel(line_wavelength)~=numel(line_weight)
    error('lsd_line_mask wavelength and weight fields have different lengths.')
end

good_lines=isfinite(line_wavelength) & isfinite(line_weight) & line_weight~=0;

line_wavelength=line_wavelength(good_lines);
line_weight=line_weight(good_lines);

if isempty(line_wavelength)
    error('No finite non-zero LSD mask lines are available.')
end

reference_weight=max(abs(line_weight));

if ~isfinite(reference_weight) || reference_weight==0
    reference_weight=1;
end

end

function [design_matrix,used_pixel_index,normalised_line_weight]=build_regularised_lsd_design_matrix(fullwave,line_wavelength,line_weight,velocity,speed_of_light,reference_weight)

%Build sparse LSD design matrix
number_wavelength_pixels=numel(fullwave);
number_velocity_points=numel(velocity);
velocity_step=median(diff(velocity),'omitnan');
velocity_min=velocity(1);
velocity_max=velocity(end);

normalised_line_weight=line_weight/reference_weight;

row_index=[];
column_index=[];
matrix_value=[];

for line_index=1:numel(line_wavelength)

    this_line_wavelength=line_wavelength(line_index);
    this_line_weight=normalised_line_weight(line_index);

    minimum_wavelength=this_line_wavelength*exp(velocity_min/speed_of_light);
    maximum_wavelength=this_line_wavelength*exp(velocity_max/speed_of_light);

    first_pixel=find(fullwave>=minimum_wavelength,1,'first');
    last_pixel=find(fullwave<=maximum_wavelength,1,'last');

    if isempty(first_pixel) || isempty(last_pixel) || last_pixel<first_pixel
        continue
    end

    pixel_index=first_pixel:last_pixel;
    relative_velocity=speed_of_light*log(fullwave(pixel_index)./this_line_wavelength);

    lower_velocity_index=floor((relative_velocity-velocity_min)/velocity_step)+1;
    upper_velocity_index=lower_velocity_index+1;

    inside_grid=lower_velocity_index>=1 & upper_velocity_index<=number_velocity_points;

    pixel_index=pixel_index(inside_grid);
    relative_velocity=relative_velocity(inside_grid);
    lower_velocity_index=lower_velocity_index(inside_grid);
    upper_velocity_index=upper_velocity_index(inside_grid);

    if isempty(pixel_index)
        continue
    end

    lower_velocity=velocity(lower_velocity_index);
    upper_weight=(relative_velocity-lower_velocity)/velocity_step;
    lower_weight=1-upper_weight;

    row_index=[row_index pixel_index pixel_index];
    column_index=[column_index lower_velocity_index upper_velocity_index];
    matrix_value=[matrix_value this_line_weight*lower_weight this_line_weight*upper_weight];

end

if isempty(row_index)
    design_matrix=sparse([],[],[],0,number_velocity_points,0);
    used_pixel_index=[];
    return
end

full_design_matrix=sparse(row_index,column_index,matrix_value,number_wavelength_pixels,number_velocity_points);
used_pixel_index=find(any(full_design_matrix~=0,2));
design_matrix=full_design_matrix(used_pixel_index,:);

end

function regularisation_matrix=make_regularisation_matrix(number_velocity_points,regularisation_order)

%Create finite-difference regularisation matrix
if regularisation_order==0
    regularisation_matrix=speye(number_velocity_points);
elseif regularisation_order==1
    e=ones(number_velocity_points,1);
    regularisation_matrix=spdiags([-e e],[0 1],number_velocity_points-1,number_velocity_points);
elseif regularisation_order==2
    e=ones(number_velocity_points,1);
    regularisation_matrix=spdiags([e -2*e e],[0 1 2],number_velocity_points-2,number_velocity_points);
else
    error('regularisation_order must be 0, 1, or 2.')
end

end

function [absorption_profile,baseline_values]=apply_lsd_baseline(velocity,absorption_profile,profile_baseline)

%Apply baseline convention to LSD absorption profiles
baseline_values=zeros(size(absorption_profile,1),1);

if strcmpi(profile_baseline,'zero')
    return
elseif strcmpi(profile_baseline,'wing_median')
    wing_region=abs(velocity)>0.8*max(abs(velocity));
    if sum(wing_region)<5
        warning('Too few points in LSD wing region. Baseline subtraction skipped.')
        return
    end
    baseline_values=median(absorption_profile(:,wing_region),2,'omitnan');
    absorption_profile=absorption_profile-repmat(baseline_values,1,size(absorption_profile,2));
else
    error(['Unknown LSD profile baseline option: ' profile_baseline])
end

end

function [absorption_profile,normalisation_scale]=normalise_lsd_absorption_profiles(velocity,absorption_profile,normalisation)

%Apply LSD normalisation convention
normalisation_scale=1;

if strcmpi(normalisation,'none')
    return
end

median_absorption=median(absorption_profile,1,'omitnan');
line_region=abs(velocity)<min(100,0.5*max(abs(velocity)));

if sum(line_region)<5
    line_region=true(size(velocity));
end

if strcmpi(normalisation,'unit_depth')
    normalisation_scale=max(abs(median_absorption(line_region)),[],'omitnan');
elseif strcmpi(normalisation,'unit_area')
    normalisation_scale=trapz(velocity,abs(median_absorption));
else
    error(['Unknown LSD normalisation option: ' normalisation])
end

if ~isfinite(normalisation_scale) || normalisation_scale==0
    warning('Could not determine a finite non-zero LSD normalisation scale. Leaving profiles unnormalised.')
    normalisation_scale=1;
end

absorption_profile=absorption_profile/normalisation_scale;

end
