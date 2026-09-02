function lsd_line_mask=make_lsd_line_mask(teff,logg,vsini,sigm,wstart,wstop,weak_limit)

%Create LSD line mask from a synthetic spectrum line list
if nargin<7 || isempty(weak_limit)
    weak_limit=1;
end

%Create synthetic line list
[synthwave,synthint,wavelengths,ele,elem,widths]=make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop);

%Force column vectors
wavelengths=wavelengths(:);
widths=widths(:);

%Validate synthetic line list
if numel(wavelengths)~=numel(widths)
    error('Synthetic line-list wavelengths and widths have different lengths.')
end

%Define retained line-mask wavelength regions
region_limits=[4260 4320;4380 4830;4895 5865;5982 6274;6326 6465;6600 6858;7415 7500];
region_labels={'S1_4260_4320','S2_4380_4830','S3_4895_5865','S4_5982_6274','S5_6326_6465','S6_6600_6858','S7_7415_7500'};

%Select lines inside retained regions
region_keep=false(size(wavelengths));
region_line_count=zeros(size(region_limits,1),1);

for region_index=1:size(region_limits,1)
    this_region=wavelengths>region_limits(region_index,1) & wavelengths<region_limits(region_index,2);
    region_keep=region_keep | this_region;
    region_line_count(region_index)=sum(this_region);
end

%Remove non-finite and weak lines
finite_keep=isfinite(wavelengths) & isfinite(widths);
strong_keep=widths>=weak_limit;
keep_lines=region_keep & finite_keep & strong_keep;

%Build selected mask arrays
selected_wavelengths=wavelengths(keep_lines);
selected_weights=widths(keep_lines);

%Sort by wavelength
[selected_wavelengths,sort_index]=sort(selected_wavelengths);
selected_weights=selected_weights(sort_index);

if isempty(selected_wavelengths)
    error('No LSD line-mask lines remain after wavelength-region and weak-line filtering.')
end

%Store mask
lsd_line_mask=struct();
lsd_line_mask.wavelength=selected_wavelengths;
lsd_line_mask.weight=selected_weights;
lsd_line_mask.width=selected_weights;

%Store synthetic spectrum metadata
lsd_line_mask.synthetic_wave=synthwave;
lsd_line_mask.synthetic_intensity=synthint;
lsd_line_mask.all_line_wavelength=wavelengths;
lsd_line_mask.all_line_weight=widths;
lsd_line_mask.element=subset_optional_line_array(ele,keep_lines,sort_index);
lsd_line_mask.element_name=subset_optional_line_array(elem,keep_lines,sort_index);

%Store selection metadata
lsd_line_mask.method='synthetic_spectrum_line_mask';
lsd_line_mask.weak_limit=weak_limit;
lsd_line_mask.region_limits=region_limits;
lsd_line_mask.region_labels=region_labels;
lsd_line_mask.region_line_count_before_weak_filter=region_line_count;
lsd_line_mask.number_lines_before_region_filter=numel(wavelengths);
lsd_line_mask.number_lines_after_region_filter=sum(region_keep & finite_keep);
lsd_line_mask.number_lines_after_weak_filter=numel(selected_wavelengths);
lsd_line_mask.teff=teff;
lsd_line_mask.logg=logg;
lsd_line_mask.vsini=vsini;
lsd_line_mask.sigm=sigm;
lsd_line_mask.wstart=wstart;
lsd_line_mask.wstop=wstop;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selected_array=subset_optional_line_array(input_array,keep_lines,sort_index)

%Subset optional line-list metadata if it matches the line-list length
selected_array=[];

if isempty(input_array)
    return
end

if numel(input_array)~=numel(keep_lines)
    return
end

if iscell(input_array)
    input_array=input_array(:);
    selected_array=input_array(keep_lines);
    selected_array=selected_array(sort_index);
elseif isstring(input_array)
    input_array=input_array(:);
    selected_array=input_array(keep_lines);
    selected_array=selected_array(sort_index);
else
    input_array=input_array(:);
    selected_array=input_array(keep_lines);
    selected_array=selected_array(sort_index);
end

end