function datacf=prepare_coarse_continuum_fit(data2,objectname,teff,logg, ...
    vsini,sigm,wstart,wstop,shift,firstorder,lastorder)
%Create a smooth synthetic-guided continuum for order merging.

output_directory=fullfile(pwd,'supplementary_data', ...
    'automatic_continuum_results');
[range_label,observed_wavelength_limits]=get_fit_range_label( ...
    data2,firstorder,lastorder);
output_file=fullfile(output_directory, ...
    [objectname '_coarse_continuum_v6_' range_label '.mat']);

if isempty(vsini), vsini=0; end

if isfile(output_file)
    loaded_fit=load(output_file,'datacf','teff','logg','vsini');
    matching_parameters=isfield(loaded_fit,'teff') && ...
        isfield(loaded_fit,'logg') && isfield(loaded_fit,'vsini') && ...
        isequaln(loaded_fit.teff,teff) && isequaln(loaded_fit.logg,logg) && ...
        isequaln(loaded_fit.vsini,vsini);
    if isfield(loaded_fit,'datacf') && matching_parameters
        fprintf('Using existing coarse automatic continuum fit\n')
        datacf=loaded_fit.datacf;
        return
    else
        fprintf('Stellar parameters changed; rebuilding coarse continuum fit\n')
    end
end

if isempty(teff) || isempty(logg)
    error('Automatic continuum fitting requires teff and logg.')
end
if ~isfolder(output_directory), mkdir(output_directory); end

config.number_anchor_bins=14;
config.synthetic_continuum_limit=0.985;
config.deep_synthetic_limit=0.75;
config.maximum_synthetic_flux=1.05;
config.minimum_points_per_bin=8;
config.edge_fraction=0.12;
config.maximum_edge_adjustment=0.05;
config.minimum_edge_points=5;
config.maximum_global_adjustment=0.10;
config.minimum_global_points=30;
config.global_continuum_percentile=99.5;
config.conservative_fallback_percentile=99.5;
config.minimum_conservative_fallback_points=20;
config.telluric_exclusion=[7590 7700];
config.telluric_anchor_width=10;
config.range_label=range_label;
config.observed_wavelength_limits=observed_wavelength_limits;

[synthetic_wave,synthetic_intensity]=make_synth_spectrum( ...
    teff,logg,vsini,sigm,wstart,wstop);
synthetic_wave=shift_velocity(synthetic_wave,shift);

datacf=struct();
fit_summary=cell(0,5);

for order_number=firstorder:lastorder
    wave_field=['wave' num2str(order_number)];
    intensity_field=['int' num2str(order_number)];
    continuum_wave_field=['w' num2str(order_number)];
    if ~isfield(data2,wave_field) || ~isfield(data2,intensity_field), continue; end

    wave=data2.(wave_field)(:).';
    representative=median(data2.(intensity_field),1,'omitnan');
    fitting_representative=representative;
    telluric_pixels=wave>=config.telluric_exclusion(1) & ...
        wave<=config.telluric_exclusion(2);
    fitting_representative(telluric_pixels)=NaN;
    synthetic=interp1(synthetic_wave,synthetic_intensity,wave, ...
        'linear',NaN);
    valid_synthetic=isfinite(synthetic) & synthetic>0 & ...
        synthetic<=config.maximum_synthetic_flux;
    synthetic_coverage=mean(valid_synthetic);

    [masked_fit,masked_info]=fit_synthetic_masked_upper_envelope( ...
        wave,fitting_representative,synthetic,config);
    [ratio_fit,ratio_info]=fit_observed_synthetic_ratio( ...
        wave,fitting_representative,synthetic,config);

    if synthetic_coverage<0.80
        [masked_fit,masked_info]=fit_observed_upper_envelope( ...
            wave,fitting_representative,config);
        ratio_info.flag=true;
    elseif masked_info.flag && ~ratio_info.flag
        masked_fit=nan(size(masked_fit));
    end

    [masked_fit,~]=stabilise_continuum_edges(wave,fitting_representative, ...
        synthetic,masked_fit,config);
    [ratio_fit,~]=stabilise_continuum_edges(wave,fitting_representative, ...
        synthetic,ratio_fit,config);
    masked_fit=stabilise_continuum_level(fitting_representative,masked_fit,config);
    ratio_fit=stabilise_continuum_level(fitting_representative,ratio_fit,config);

    [continuum,method]=combine_candidate_fits(masked_fit,ratio_fit, ...
        masked_info,ratio_info);
    continuum=stabilise_continuum_level(fitting_representative,continuum,config);
    valid_continuum=isfinite(continuum) & continuum>0;

    if mean(valid_continuum)<0.80
        [continuum,~]=fit_observed_upper_envelope( ...
            wave,fitting_representative,config);
        continuum=stabilise_continuum_level( ...
            fitting_representative,continuum,config);
        method='observed_upper_envelope';
        valid_continuum=isfinite(continuum) & continuum>0;
    end
    if mean(valid_continuum)<0.80
        continuum=fit_conservative_constant_upper_envelope( ...
            fitting_representative,config);
        method='constant_upper_envelope_fallback';
        valid_continuum=isfinite(continuum) & continuum>0;
        if mean(valid_continuum)>=0.80
            fprintf(['Using conservative constant upper-envelope ' ...
                'fallback for order %.0f\n'],order_number)
        end
    end
    if mean(valid_continuum)<0.80
        continuum=nan(size(wave));
        method='unusable_order_zero_weight';
        valid_continuum=false(size(wave));
        fprintf(['Order %.0f has insufficient positive data for a ' ...
            'continuum; excluding it from the merge\n'],order_number)
    end
    continuum=bridge_telluric_continuum(wave,continuum,config);

    datacf.(continuum_wave_field)=wave;
    datacf.(intensity_field)=continuum;
    fit_summary(end+1,:)={order_number,method,synthetic_coverage, ...
        mean(valid_continuum),median(representative./continuum,'omitnan')}; %#ok<AGROW>
end

datacf.automatic_fit_method=['synthetic_guided_hybrid_v6_' ...
    range_label '_o2_masked'];
datacf.automatic_fit_summary=cell2table(fit_summary,'VariableNames', ...
    {'Order','Method','SyntheticCoverage','ContinuumCoverage', ...
    'MedianNormalisedIntensity'});
save(output_file,'datacf','config','teff','logg','vsini','-v7.3')

end

function continuum=fit_conservative_constant_upper_envelope(flux,config)

%Last-resort Stage 1 fit that cannot follow broad absorption structure
continuum=nan(size(flux));
valid_flux=isfinite(flux) & flux>0;
if sum(valid_flux)<config.minimum_conservative_fallback_points, return; end

upper_level=prctile(flux(valid_flux), ...
    config.conservative_fallback_percentile);
if ~isfinite(upper_level) || upper_level<=0, return; end
continuum(:)=upper_level;

end

function [range_label,wavelength_limits]=get_fit_range_label( ...
    data2,firstorder,lastorder)

minimum_wavelength=Inf;
maximum_wavelength=-Inf;
for order_number=firstorder:lastorder
    wave_field=['wave' num2str(order_number)];
    if ~isfield(data2,wave_field) || isempty(data2.(wave_field)), continue; end
    minimum_wavelength=min(minimum_wavelength,min(data2.(wave_field)));
    maximum_wavelength=max(maximum_wavelength,max(data2.(wave_field)));
end

if ~isfinite(minimum_wavelength) || ~isfinite(maximum_wavelength)
    error('No wavelength data are available for coarse continuum fitting.')
end

wavelength_limits=[minimum_wavelength maximum_wavelength];
if minimum_wavelength<3860 || maximum_wavelength>6875
    range_label='extended';
else
    range_label='standard';
end

end

function continuum=bridge_telluric_continuum(wave,continuum,config)

%Do not let the coarse continuum follow the saturated O2 A-band
region=config.telluric_exclusion;
inside=wave>=region(1) & wave<=region(2);
if ~any(inside), return; end

left=wave>=region(1)-config.telluric_anchor_width & ...
    wave<region(1) & isfinite(continuum) & continuum>0;
right=wave>region(2) & ...
    wave<=region(2)+config.telluric_anchor_width & ...
    isfinite(continuum) & continuum>0;

if any(left) && any(right)
    anchor_values=[median(continuum(left),'omitnan') ...
        median(continuum(right),'omitnan')];
    continuum(inside)=interp1(region,anchor_values,wave(inside),'linear');
elseif any(left)
    continuum(inside)=median(continuum(left),'omitnan');
elseif any(right)
    continuum(inside)=median(continuum(right),'omitnan');
else
    continuum(inside)=NaN;
end

end

function [continuum,info]=fit_synthetic_masked_upper_envelope( ...
    wave,flux,synthetic,config)

synthetic_gradient=abs(gradient(synthetic)./max(abs(gradient(wave)),eps));
finite_gradient=synthetic_gradient(isfinite(synthetic_gradient));
if isempty(finite_gradient)
    gradient_limit=Inf;
else
    gradient_limit=prctile(finite_gradient,70);
end
candidate=isfinite(wave) & isfinite(flux) & flux>0 & ...
    isfinite(synthetic) & synthetic>=config.synthetic_continuum_limit & ...
    synthetic<=config.maximum_synthetic_flux & ...
    synthetic_gradient<=gradient_limit;
[anchor_wave,anchor_flux]=make_upper_envelope_anchors(wave,flux,candidate, ...
    config.number_anchor_bins,config.minimum_points_per_bin);
[continuum,anchor_wave,anchor_flux]=iterative_anchor_fit( ...
    wave,flux,anchor_wave,anchor_flux,candidate);
info=make_fit_info(continuum,anchor_wave,anchor_flux);

end

function [continuum,info]=fit_observed_synthetic_ratio( ...
    wave,flux,synthetic,config)

valid=isfinite(wave) & isfinite(flux) & flux>0 & isfinite(synthetic) & ...
    synthetic>config.deep_synthetic_limit & ...
    synthetic<=config.maximum_synthetic_flux;
ratio=nan(size(flux));
ratio(valid)=flux(valid)./synthetic(valid);
[anchor_wave,anchor_ratio]=make_robust_median_anchors(wave,ratio,valid, ...
    config.number_anchor_bins,config.minimum_points_per_bin);
continuum=fit_smooth_positive_polynomial(wave,anchor_wave,anchor_ratio);
info=make_fit_info(continuum,anchor_wave,anchor_ratio);

end

function [continuum,info]=fit_observed_upper_envelope(wave,flux,config)

candidate=isfinite(wave) & isfinite(flux) & flux>0;
[anchor_wave,anchor_flux]=make_upper_envelope_anchors(wave,flux,candidate, ...
    config.number_anchor_bins,config.minimum_points_per_bin);
[fit_wave,fit_flux]=select_upper_hull_anchors(anchor_wave,anchor_flux);
if numel(fit_wave)<4
    fit_wave=anchor_wave;
    fit_flux=anchor_flux;
end
continuum=fit_smooth_positive_polynomial(wave,fit_wave,fit_flux);
info=make_fit_info(continuum,anchor_wave,anchor_flux);

end

function [continuum,method]=combine_candidate_fits( ...
    masked_fit,ratio_fit,masked_info,ratio_info)

valid_masked=isfinite(masked_fit) & masked_fit>0;
valid_ratio=isfinite(ratio_fit) & ratio_fit>0;
masked_coverage=mean(valid_masked);
ratio_coverage=mean(valid_ratio);
continuum=nan(size(masked_fit));
both=valid_masked & valid_ratio;

if masked_coverage>=0.80 && ratio_coverage>=0.80 && ...
        ~masked_info.flag && ~ratio_info.flag
    continuum(both)=sqrt(masked_fit(both).*ratio_fit(both));
    continuum(valid_masked & ~valid_ratio)=masked_fit(valid_masked & ~valid_ratio);
    continuum(valid_ratio & ~valid_masked)=ratio_fit(valid_ratio & ~valid_masked);
    method='hybrid';
elseif ratio_coverage>=0.80 && ~ratio_info.flag
    continuum=ratio_fit;
    method='synthetic_ratio';
elseif masked_coverage>=0.80
    continuum=masked_fit;
    method='synthetic_masked';
elseif ratio_coverage>=0.80
    continuum=ratio_fit;
    method='synthetic_ratio_review';
else
    method='failed';
end

end

function [anchor_wave,anchor_flux]=make_upper_envelope_anchors( ...
    wave,flux,candidate,number_bins,minimum_points)

edges=linspace(min(wave),max(wave),number_bins+1);
anchor_wave=[];
anchor_flux=[];
for bin_index=1:number_bins
    in_bin=candidate & wave>=edges(bin_index) & wave<=edges(bin_index+1);
    if nnz(in_bin)<minimum_points, continue; end
    upper_limit=prctile(flux(in_bin),80);
    upper=in_bin & flux>=upper_limit;
    anchor_wave(end+1)=median(wave(upper),'omitnan'); %#ok<AGROW>
    anchor_flux(end+1)=median(flux(upper),'omitnan'); %#ok<AGROW>
end

end

function [anchor_wave,anchor_value]=make_robust_median_anchors( ...
    wave,value,valid,number_bins,minimum_points)

edges=linspace(min(wave),max(wave),number_bins+1);
anchor_wave=[];
anchor_value=[];
for bin_index=1:number_bins
    in_bin=valid & wave>=edges(bin_index) & wave<=edges(bin_index+1);
    if nnz(in_bin)<minimum_points, continue; end
    limits=prctile(value(in_bin),[20 80]);
    robust=in_bin & value>=limits(1) & value<=limits(2);
    anchor_wave(end+1)=median(wave(robust),'omitnan'); %#ok<AGROW>
    anchor_value(end+1)=median(value(robust),'omitnan'); %#ok<AGROW>
end

end

function [hull_wave,hull_flux]=select_upper_hull_anchors( ...
    anchor_wave,anchor_flux)

valid=isfinite(anchor_wave) & isfinite(anchor_flux);
[anchor_wave,sort_index]=sort(anchor_wave(valid));
anchor_flux=anchor_flux(valid);
anchor_flux=anchor_flux(sort_index);
hull_index=[];
for anchor_index=1:numel(anchor_wave)
    hull_index(end+1)=anchor_index; %#ok<AGROW>
    while numel(hull_index)>=3
        first=hull_index(end-2);
        middle=hull_index(end-1);
        last=hull_index(end);
        first_slope=(anchor_flux(middle)-anchor_flux(first))/ ...
            max(anchor_wave(middle)-anchor_wave(first),eps);
        second_slope=(anchor_flux(last)-anchor_flux(middle))/ ...
            max(anchor_wave(last)-anchor_wave(middle),eps);
        if second_slope>first_slope
            hull_index(end-1)=[];
        else
            break
        end
    end
end
hull_wave=anchor_wave(hull_index);
hull_flux=anchor_flux(hull_index);

end

function [continuum,anchor_wave,anchor_flux]=iterative_anchor_fit( ...
    wave,flux,anchor_wave,anchor_flux,candidate)

if numel(anchor_wave)<4
    continuum=nan(size(wave));
    return
end
continuum=pchip_with_endpoint_hold(anchor_wave,anchor_flux,wave);
for iteration=1:5
    residual=flux./continuum-1;
    scale=robust_sigma(residual(candidate));
    if ~isfinite(scale) || scale<=0, break; end
    refined=candidate & residual>-1.5*scale & residual<3*scale;
    [new_wave,new_flux]=make_upper_envelope_anchors(wave,flux,refined, ...
        max(8,numel(anchor_wave)),5);
    if numel(new_wave)<4, break; end
    new_continuum=pchip_with_endpoint_hold(new_wave,new_flux,wave);
    change=median(abs(new_continuum./continuum-1),'omitnan');
    continuum=new_continuum;
    anchor_wave=new_wave;
    anchor_flux=new_flux;
    if change<1e-4, break; end
end
continuum(~isfinite(continuum) | continuum<=0)=NaN;

end

function continuum=fit_smooth_positive_polynomial( ...
    wave,anchor_wave,anchor_value)

valid=isfinite(anchor_wave) & isfinite(anchor_value) & anchor_value>0;
anchor_wave=anchor_wave(valid);
anchor_value=anchor_value(valid);
if numel(anchor_wave)<4
    continuum=nan(size(wave));
    return
end
keep=true(size(anchor_wave));
for iteration=1:5
    polynomial_order=min(3,nnz(keep)-1);
    [polynomial,fit_structure,centre_scale]=polyfit(anchor_wave(keep), ...
        log(anchor_value(keep)),polynomial_order);
    fitted=exp(polyval(polynomial,anchor_wave,fit_structure,centre_scale));
    residual=log(anchor_value./fitted);
    scale=robust_sigma(residual(keep));
    if ~isfinite(scale) || scale<=0, break; end
    new_keep=abs(residual)<=3*scale;
    if nnz(new_keep)<4 || isequal(new_keep,keep), break; end
    keep=new_keep;
end
polynomial_order=min(3,nnz(keep)-1);
[polynomial,fit_structure,centre_scale]=polyfit(anchor_wave(keep), ...
    log(anchor_value(keep)),polynomial_order);
continuum=exp(polyval(polynomial,wave,fit_structure,centre_scale));
support_min=min(anchor_wave(keep));
support_max=max(anchor_wave(keep));
continuum(wave<support_min)=exp(polyval(polynomial,support_min, ...
    fit_structure,centre_scale));
continuum(wave>support_max)=exp(polyval(polynomial,support_max, ...
    fit_structure,centre_scale));
continuum(~isfinite(continuum) | continuum<=0)=NaN;

end

function interpolated=pchip_with_endpoint_hold( ...
    input_wave,input_value,output_wave)

valid=isfinite(input_wave) & isfinite(input_value) & input_value>0;
input_wave=input_wave(valid);
input_value=input_value(valid);
[input_wave,sort_index]=sort(input_wave);
input_value=input_value(sort_index);
[input_wave,unique_index]=unique(input_wave,'stable');
input_value=input_value(unique_index);
if numel(input_wave)<2
    interpolated=nan(size(output_wave));
    return
end
interpolated=interp1(input_wave,input_value,output_wave,'pchip',NaN);
interpolated(output_wave<input_wave(1))=input_value(1);
interpolated(output_wave>input_wave(end))=input_value(end);

end

function [continuum,edge_info]=stabilise_continuum_edges( ...
    wave,flux,synthetic,continuum,config)

edge_info=struct('blue_scale',NaN,'red_scale',NaN);
valid_continuum=isfinite(continuum) & continuum>0;
if mean(valid_continuum)<0.80, return; end
trusted=isfinite(wave) & isfinite(flux) & flux>0 & ...
    isfinite(synthetic) & synthetic>=config.synthetic_continuum_limit & ...
    synthetic<=config.maximum_synthetic_flux & valid_continuum;
wave_span=max(wave)-min(wave);
blue_limit=min(wave)+config.edge_fraction*wave_span;
red_limit=max(wave)-config.edge_fraction*wave_span;
blue=trusted & wave<=blue_limit;
red=trusted & wave>=red_limit;
blue_scale=1;
red_scale=1;
if nnz(blue)>=config.minimum_edge_points
    blue_scale=median(flux(blue)./continuum(blue),'omitnan');
end
if nnz(red)>=config.minimum_edge_points
    red_scale=median(flux(red)./continuum(red),'omitnan');
end
minimum_scale=1-config.maximum_edge_adjustment;
maximum_scale=1+config.maximum_edge_adjustment;
blue_scale=min(max(blue_scale,minimum_scale),maximum_scale);
red_scale=min(max(red_scale,minimum_scale),maximum_scale);
edge_info.blue_scale=blue_scale;
edge_info.red_scale=red_scale;

edge_correction=ones(size(wave));
blue_position=(wave-min(wave))./max(blue_limit-min(wave),eps);
blue_position=min(max(blue_position,0),1);
blue_smooth=blue_position.^2.*(3-2*blue_position);
blue_region=wave<=blue_limit;
edge_correction(blue_region)=blue_scale+(1-blue_scale).* ...
    blue_smooth(blue_region);
red_position=(max(wave)-wave)./max(max(wave)-red_limit,eps);
red_position=min(max(red_position,0),1);
red_smooth=red_position.^2.*(3-2*red_position);
red_region=wave>=red_limit;
edge_correction(red_region)=red_scale+(1-red_scale).* ...
    red_smooth(red_region);
continuum=continuum.*edge_correction;

end

function continuum=stabilise_continuum_level(flux,continuum,config)

trusted=isfinite(flux) & flux>0 & isfinite(continuum) & continuum>0;
if nnz(trusted)<config.minimum_global_points, return; end
scale=prctile(flux(trusted)./continuum(trusted), ...
    config.global_continuum_percentile);
scale=min(max(scale,1-config.maximum_global_adjustment), ...
    1+config.maximum_global_adjustment);
continuum=continuum.*scale;

end

function info=make_fit_info(continuum,anchor_wave,anchor_flux)

info.anchor_wave=anchor_wave;
info.anchor_flux=anchor_flux;
info.number_anchors=numel(anchor_wave);
valid=isfinite(continuum) & continuum>0;
if nnz(valid)>=5
    excessive_curvature=max(abs(diff(log(continuum(valid)),2)), ...
        [],'omitnan')>0.02;
else
    excessive_curvature=true;
end
info.flag=mean(valid)<0.95 || ...
    (info.number_anchors>0 && info.number_anchors<6) || excessive_curvature;

end

function sigma=robust_sigma(values)

values=values(isfinite(values));
if isempty(values)
    sigma=NaN;
else
    centre=median(values);
    sigma=1.4826*median(abs(values-centre));
end

end
