function [fullint,cleanup]=remove_isolated_merged_cosmic_ray_spikes( ...
    fullwave,fullint)
%REMOVE_ISOLATED_MERGED_COSMIC_RAY_SPIKES Light-touch post-merge cleanup.
% A candidate must be a narrow positive feature that is significant both
% relative to its local spectrum and relative to the other observations.

local_baseline_window=31;
local_noise_window=101;
local_sigma_threshold=8;
temporal_sigma_threshold=6;
support_sigma_threshold=4;
maximum_event_width=3;
edge_trim=31;

cleanup=struct();
cleanup.method='light_touch_positive_local_and_temporal_spike';
cleanup.local_baseline_window=local_baseline_window;
cleanup.local_noise_window=local_noise_window;
cleanup.local_sigma_threshold=local_sigma_threshold;
cleanup.temporal_sigma_threshold=temporal_sigma_threshold;
cleanup.support_sigma_threshold=support_sigma_threshold;
cleanup.maximum_event_width=maximum_event_width;
cleanup.edge_trim=edge_trim;
cleanup.number_pixels=0;
cleanup.number_events=0;
cleanup.observation_index=[];
cleanup.wavelength=[];
cleanup.original_intensity=[];
cleanup.replacement_intensity=[];

if size(fullint,1)<3
    cleanup.status='skipped: fewer than three spectra';
    fprintf(['Skipping light-touch merged-spectrum cosmic-ray filtering: ' ...
        'at least three spectra are required\n'])
    return
end

if size(fullint,2)<=2*edge_trim || numel(fullwave)~=size(fullint,2)
    cleanup.status='skipped: incompatible or short wavelength grid';
    warning(['Skipping light-touch merged-spectrum cosmic-ray filtering: ' ...
        'the wavelength grid is incompatible or too short.'])
    return
end

local_baseline=movmedian(fullint,local_baseline_window,2,'omitnan');
local_residual=fullint-local_baseline;
local_centre=movmedian(local_residual,local_noise_window,2,'omitnan');
local_sigma=1.4826*movmedian(abs(local_residual-local_centre), ...
    local_noise_window,2,'omitnan');

row_noise_floor=0.25*median(local_sigma,2,'omitnan');
row_fallback=std(local_residual,0,2,'omitnan');
invalid_floor=~isfinite(row_noise_floor) | row_noise_floor<=0;
row_noise_floor(invalid_floor)=row_fallback(invalid_floor);
row_noise_floor(~isfinite(row_noise_floor) | row_noise_floor<=0)=eps;
noise_floor_image=repmat(row_noise_floor,1,size(fullint,2));
invalid_sigma=~isfinite(local_sigma) | local_sigma<=0;
local_sigma(invalid_sigma)=noise_floor_image(invalid_sigma);
local_sigma=max(local_sigma,noise_floor_image);
local_significance=local_residual./local_sigma;

temporal_centre=median(fullint,1,'omitnan');
temporal_residual=fullint-temporal_centre;
temporal_sigma=1.4826*median(abs(temporal_residual),1,'omitnan');
temporal_noise_floor=0.5*median(local_sigma,1,'omitnan');
temporal_noise_floor(~isfinite(temporal_noise_floor) | ...
    temporal_noise_floor<=0)=eps;
invalid_temporal_sigma=~isfinite(temporal_sigma) | temporal_sigma<=0;
temporal_sigma(invalid_temporal_sigma)= ...
    temporal_noise_floor(invalid_temporal_sigma);
temporal_sigma=max(temporal_sigma,temporal_noise_floor);

valid_pixel=isfinite(fullint) & isfinite(local_baseline) & ...
    isfinite(local_significance);
candidate_seed=local_significance>local_sigma_threshold & ...
    temporal_residual>temporal_sigma_threshold*temporal_sigma & ...
    valid_pixel;
candidate_support=local_significance>support_sigma_threshold & ...
    local_residual>0 & valid_pixel;
candidate_seed(:,1:edge_trim)=false;
candidate_seed(:,end-edge_trim+1:end)=false;
candidate_support(:,1:edge_trim)=false;
candidate_support(:,end-edge_trim+1:end)=false;
candidate_mask=grow_narrow_events(candidate_seed,candidate_support, ...
    maximum_event_width);

original_intensity=fullint;
fullint(candidate_mask)=local_baseline(candidate_mask);
[observation_index,pixel_index]=find(candidate_mask);
cleanup.status='completed';
cleanup.number_pixels=numel(observation_index);
cleanup.number_events=count_events(candidate_mask);
cleanup.observation_index=observation_index;
cleanup.wavelength=fullwave(pixel_index(:));
cleanup.original_intensity=original_intensity(candidate_mask);
cleanup.replacement_intensity=fullint(candidate_mask);

fprintf(['Light-touch merged-spectrum cosmic-ray filter replaced %.0f ' ...
    'pixels in %.0f narrow events\n'],cleanup.number_pixels, ...
    cleanup.number_events)

end

function output_mask=grow_narrow_events(seed,support,maximum_width)

output_mask=false(size(seed));
for row_index=1:size(seed,1)
    padded=[false support(row_index,:) false];
    run_starts=find(diff(padded)==1);
    run_stops=find(diff(padded)==-1)-1;
    for run_index=1:numel(run_starts)
        indices=run_starts(run_index):run_stops(run_index);
        if numel(indices)<=maximum_width && any(seed(row_index,indices))
            output_mask(row_index,indices)=true;
        end
    end
end

end

function number_events=count_events(mask)

number_events=0;
for row_index=1:size(mask,1)
    number_events=number_events+sum(diff([false mask(row_index,:) false])==1);
end

end
