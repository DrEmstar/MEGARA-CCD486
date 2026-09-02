function [normalised_intensity,stage2]=apply_suppnet_merged_continuum( ...
    wavelength,intensity,objectname)
%Apply one SUPPNet continuum correction to a complete merged spectrum.

output_directory=fullfile(pwd,'supplementary_data', ...
    'automatic_continuum_results');
if ~isfolder(output_directory), mkdir(output_directory); end

wavelength=wavelength(:).';
if size(intensity,2)~=numel(wavelength)
    error('Merged wavelength and intensity dimensions do not match.')
end

representative=median(intensity,1,'omitnan');
telluric_region=[7590 7700];
[suppnet_representative,telluric_bridge]=bridge_masked_region( ...
    wavelength,representative,telluric_region);
valid_input=isfinite(wavelength) & isfinite(suppnet_representative) & ...
    suppnet_representative>0;
if mean(valid_input)<0.95
    warning('Only %.1f%% of the merged spectrum is valid for SUPPNet.', ...
        100*mean(valid_input))
end
if nnz(valid_input)<1000
    error('The merged spectrum contains too few valid pixels for SUPPNet.')
end

input_file=fullfile(output_directory, ...
    [objectname '_merged_stage2.txt']);
writematrix([wavelength(valid_input).' suppnet_representative(valid_input).'], ...
    input_file,'Delimiter',' ','FileType','text')

suppnet_options=struct('quiet',true,'sampling',[], ...
    'weights','synth','skip',0,'extraArgs','');
run_suppnet_spectrum(input_file,suppnet_options)
suppnet_result=load_suppnet_outputs(input_file);

suppnet_wave=suppnet_result.wave(:).';
suppnet_fit=suppnet_result.fit(:).';
valid_suppnet=isfinite(suppnet_wave) & isfinite(suppnet_fit) & ...
    suppnet_fit>0;
suppnet_wave=suppnet_wave(valid_suppnet);
suppnet_fit=suppnet_fit(valid_suppnet);

input_median=median(suppnet_representative(valid_input),'omitnan');
continuum=input_median*interp1(suppnet_wave,suppnet_fit, ...
    wavelength,'linear',NaN);
[continuum,continuum_bridge]=bridge_masked_region( ...
    wavelength,continuum,telluric_region);
valid_continuum=isfinite(continuum) & continuum>0;
if mean(valid_continuum(valid_input))<0.95
    error('SUPPNet continuum covers less than 95%% of the valid merged spectrum.')
end

continuum_limits=prctile(continuum(valid_continuum),[0.1 99.9]);
if continuum_limits(1)<0.25 || continuum_limits(2)>4
    error(['SUPPNet returned an implausible merged continuum range: ' ...
        num2str(continuum_limits(1)) ' to ' num2str(continuum_limits(2)) '.'])
end

normalised_intensity=intensity./continuum;
invalid=repmat(~valid_continuum,size(intensity,1),1) | ...
    ~isfinite(intensity);
normalised_intensity(invalid)=NaN;

stage2.method='SUPPNet merged-spectrum correction';
stage2.input_file=input_file;
stage2.output_file=suppnet_result.outputFile;
stage2.wavelength=wavelength;
stage2.representative_input=representative;
stage2.suppnet_representative=suppnet_representative;
stage2.continuum=continuum;
stage2.valid_continuum=valid_continuum;
stage2.input_median=input_median;
stage2.valid_coverage=mean(valid_continuum(valid_input));
stage2.continuum_limits=continuum_limits;
stage2.telluric_bridge_region=telluric_region;
stage2.telluric_input_bridge=telluric_bridge;
stage2.telluric_continuum_bridge=continuum_bridge;

save(fullfile(output_directory, ...
    [objectname '_merged_stage2_continuum.mat']),'stage2','-v7.3')

end

function [values,bridge]=bridge_masked_region(wavelength,values,region)

%Replace the O2 band with a smooth bridge between nearby continuum anchors
anchor_width=10;
inside=wavelength>=region(1) & wavelength<=region(2);
if ~any(inside)
    bridge.left_value=NaN;
    bridge.right_value=NaN;
    bridge.anchor_width=anchor_width;
    bridge.number_replaced=0;
    return
end

left_anchor=wavelength>=region(1)-anchor_width & ...
    wavelength<region(1) & isfinite(values);
right_anchor=wavelength>region(2) & ...
    wavelength<=region(2)+anchor_width & isfinite(values);

if nnz(left_anchor)<5 || nnz(right_anchor)<5
    error('Insufficient valid data to bridge the O2 A-band continuum.')
end

left_value=median(values(left_anchor),'omitnan');
right_value=median(values(right_anchor),'omitnan');
values(inside)=interp1(region,[left_value right_value], ...
    wavelength(inside),'linear');

bridge.left_value=left_value;
bridge.right_value=right_value;
bridge.anchor_width=anchor_width;
bridge.number_replaced=sum(inside);

end
