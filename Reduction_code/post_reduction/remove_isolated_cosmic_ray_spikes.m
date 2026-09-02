function [data2,cleanup_report]=remove_isolated_cosmic_ray_spikes( ...
    data2,firstorder,lastorder)
%REMOVE_ISOLATED_COSMIC_RAY_SPIKES Clean conservative 1-D cosmic candidates.
% Candidates must be positive, locally significant, temporally unusual
% when at least three spectra are present, and narrow.  Noise is estimated
% locally along each order so a noisy part of an order does not make the
% filter insensitive everywhere else.

fprintf('Removing isolated positive cosmic-ray spikes\n')

local_window=31;
local_noise_window=101;
edge_trim=31;
local_sigma_threshold=6;
support_sigma_threshold=3;
temporal_sigma_threshold=5;
extreme_local_sigma_threshold=12;
maximum_run_width=5;

cleanup_report=struct();
cleanup_report.method='positive_local_and_temporal_spike';
cleanup_report.local_window=local_window;
cleanup_report.local_noise_window=local_noise_window;
cleanup_report.local_sigma_threshold=local_sigma_threshold;
cleanup_report.support_sigma_threshold=support_sigma_threshold;
cleanup_report.temporal_sigma_threshold=temporal_sigma_threshold;
cleanup_report.extreme_local_sigma_threshold= ...
    extreme_local_sigma_threshold;
cleanup_report.maximum_run_width=maximum_run_width;
cleanup_report.edge_trim=edge_trim;
cleanup_report.orders=struct('order',{},'candidate_count',{}, ...
    'observation_index',{},'wavelength',{},'original_intensity',{}, ...
    'replacement_intensity',{});
cleanup_report.total_candidate_count=0;

report_index=0;

for ord=firstorder:lastorder

    wave_field_name=['wave' num2str(ord)];
    intensity_field_name=['int' num2str(ord)];
    weights_field_name=['weights' num2str(ord)];

    if ~isfield(data2,intensity_field_name) || ...
            ~isfield(data2,wave_field_name) || ...
            ~isfield(data2,weights_field_name)
        continue
    end

    int=data2.(intensity_field_name);
    wave=data2.(wave_field_name);
    weights=data2.(weights_field_name);

    if size(int,2)<=2*edge_trim
        warning(['Order ' num2str(ord) ' is too short for cosmic-ray ' ...
            'filtering and edge trimming.'])
        continue
    end

    local_baseline=movmedian(int,local_window,2,'omitnan');
    local_residual=int-local_baseline;
    local_centre=movmedian(local_residual,local_noise_window,2, ...
        'omitnan');
    local_sigma=1.4826*movmedian(abs(local_residual-local_centre), ...
        local_noise_window,2,'omitnan');

    %Prevent an exceptionally smooth local interval from producing an
    %unrealistically small threshold.  The floor is calculated separately
    %for each observation and is deliberately only a fraction of its
    %typical rolling noise.
    row_noise_floor=0.25*median(local_sigma,2,'omitnan');
    local_fallback=std(local_residual,0,2,'omitnan');
    invalid_row_floor=~isfinite(row_noise_floor) | row_noise_floor<=0;
    row_noise_floor(invalid_row_floor)=local_fallback(invalid_row_floor);
    row_noise_floor(~isfinite(row_noise_floor) | row_noise_floor<=0)=eps;
    invalid_local_sigma=~isfinite(local_sigma) | local_sigma<=0;
    row_noise_floor_image=repmat(row_noise_floor,1,size(int,2));
    local_sigma(invalid_local_sigma)=row_noise_floor_image( ...
        invalid_local_sigma);
    local_sigma=max(local_sigma,row_noise_floor_image);
    local_significance=local_residual./local_sigma;
    local_candidate=local_significance>local_sigma_threshold;

    if size(int,1)>=3
        temporal_centre=median(int,1,'omitnan');
        temporal_residual=int-temporal_centre;
        temporal_sigma=1.4826*median( ...
            abs(temporal_residual),1,'omitnan');
        invalid_temporal_sigma=~isfinite(temporal_sigma) | ...
            temporal_sigma<=0;
        temporal_noise_floor=0.5*median(local_sigma,1,'omitnan');
        temporal_noise_floor(~isfinite(temporal_noise_floor) | ...
            temporal_noise_floor<=0)=eps;
        temporal_sigma(invalid_temporal_sigma)= ...
            temporal_noise_floor(invalid_temporal_sigma);
        temporal_sigma=max(temporal_sigma,temporal_noise_floor);
        temporal_sigma(~isfinite(temporal_sigma) | ...
            temporal_sigma<=0)=Inf;
        temporal_candidate=temporal_residual > ...
            temporal_sigma_threshold*temporal_sigma;
        candidate_seed=local_candidate & temporal_candidate;
    else
        candidate_seed=local_significance> ...
            extreme_local_sigma_threshold;
    end

    valid_pixel=isfinite(int) & isfinite(local_baseline) & ...
        isfinite(local_significance);
    candidate_seed=candidate_seed & valid_pixel;
    candidate_support=local_residual>0 & ...
        local_significance>support_sigma_threshold & valid_pixel;
    candidate_seed(:,1:edge_trim)=false;
    candidate_seed(:,end-edge_trim+1:end)=false;
    candidate_support(:,1:edge_trim)=false;
    candidate_support(:,end-edge_trim+1:end)=false;
    candidate_mask=grow_narrow_positive_events(candidate_seed, ...
        candidate_support,maximum_run_width);

    original_intensity=int;
    int(candidate_mask)=local_baseline(candidate_mask);

    report_index=report_index+1;
    [observation_index,pixel_index]=find(candidate_mask);
    cleanup_report.orders(report_index).order=ord;
    cleanup_report.orders(report_index).candidate_count= ...
        numel(observation_index);
    cleanup_report.orders(report_index).observation_index= ...
        observation_index;
    cleanup_report.orders(report_index).wavelength=wave(pixel_index(:));
    cleanup_report.orders(report_index).original_intensity= ...
        original_intensity(candidate_mask);
    cleanup_report.orders(report_index).replacement_intensity= ...
        int(candidate_mask);
    cleanup_report.total_candidate_count= ...
        cleanup_report.total_candidate_count+numel(observation_index);

    data2.(intensity_field_name)=int(:,edge_trim+1:end-edge_trim);
    data2.(wave_field_name)=wave(edge_trim+1:end-edge_trim);
    data2.(weights_field_name)=weights(:,edge_trim+1:end-edge_trim);

end


fprintf('Replaced %.0f isolated positive cosmic-ray candidates\n', ...
    cleanup_report.total_candidate_count)

end


function output_mask=grow_narrow_positive_events( ...
    candidate_seed,candidate_support,maximum_run_width)

output_mask=false(size(candidate_seed));

for row_index=1:size(candidate_seed,1)

    padded_mask=[false candidate_support(row_index,:) false];
    run_starts=find(diff(padded_mask)==1);
    run_stops=find(diff(padded_mask)==-1)-1;

    for run_index=1:numel(run_starts)

        run_width=run_stops(run_index)-run_starts(run_index)+1;

        run_contains_seed=any(candidate_seed(row_index, ...
            run_starts(run_index):run_stops(run_index)));

        if run_contains_seed && run_width<=maximum_run_width
            output_mask(row_index,run_starts(run_index): ...
                run_stops(run_index))=true;
        end

    end

end


end
