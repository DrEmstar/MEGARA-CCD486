function [myflag,finaldata]=post_reduction(objectname,options,teff,logg,vsini)

%DEFINITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigm=0.5;wstart=3750;wstop=8950;
%systemic velocity
shift=0;
firstorder=55;
lastorder=154;
step=0.02;

%Initialise outputs
myflag=1;
finaldata=[];

input_vsini=vsini;
synthetic_vsini=get_synthetic_vsini(input_vsini);

%Initialise fallback measurements
radial_velocity=NaN;
radial_velocity_error=NaN;
vsini=NaN;
vsini_error=NaN;
signal_to_noise=NaN;

target_directory=pwd;
reduced_frames_directory=fullfile(target_directory,'reduced frames');
processing_directory=fullfile(reduced_frames_directory,'processing_files');
supplementary_directory=fullfile(target_directory,'supplementary_data');
if ~isfolder(reduced_frames_directory), mkdir(reduced_frames_directory); end
if ~isfolder(processing_directory), mkdir(processing_directory); end
if ~isfolder(supplementary_directory), mkdir(supplementary_directory); end

final_data_file_name=get_final_data_file_name(objectname,options);

weak_limit1=1; %lines < XX mA in strength removed from LSD line mask (default 1)

if ~isfield(options,'calculate_lsd_profile')
    error('Missing options.calculate_lsd_profile. MEGARA 2.1 requires this option for regularised LSD profile control.')
end

if ~isfield(options,'chunked_post_reduction')
    options.chunked_post_reduction=false;
end

if ~isfield(options,'post_reduction_chunk_size') || isempty(options.post_reduction_chunk_size)
    options.post_reduction_chunk_size=500;
end
if ~isfield(options,'extended_wavelength_range') || ...
        isempty(options.extended_wavelength_range)
    options.extended_wavelength_range=false;
end

if options.extended_wavelength_range
    warning(['Extended wavelength range enabled: approximately 3750-8950 ' ...
        'Angstrom. The extreme blue and red regions may require manual ' ...
        'review and additional processing.'])
else
    fprintf('Using the standard wavelength range: 3869-6866 Angstrom\n')
end

%list files
list=list_reduced_spectrum_files(reduced_frames_directory);
if isempty(list), warning(['No reduced J*.mat files found for ' objectname '.']); myflag=1; finaldata=[]; return; end
list=sort_reduced_spectrum_files(list);
%Load first spectrum
test_data=load_primary_variable_from_mat_file(get_list_file_path(list(1)));

%checks maximum order that has been reduced
g=lastorder;
order_exist_test=false;
while ~order_exist_test && g>=firstorder
    order_exist_test=isfield(test_data,['order_' num2str(g)]);
    g=g-1;
end
if ~order_exist_test, error(['No order fields found between order_' num2str(firstorder) ' and order_' num2str(lastorder) ' in ' list(1).name '.']); end
lastorder=g+1;

%Check for previous final data
test1=dir(final_data_file_name);

%Create or load full data
if numel(test1)==0 || (numel(test1)>0 && options.overwrite_full_data==1)

    if options.apply_continuum_fit~=1
        error('Cannot create final_data without apply_continuum_fit unless a previous final_data file already exists.')
    end

    if options.order_merge~=1
        error('Chunked post-reduction currently requires options.order_merge=true.')
    end

    if options.manual_continuum_fit==1 && options.chunked_post_reduction==1
        error('Manual continuum fitting cannot be created during chunked post-reduction. Create datacf first, then rerun with options.manual_continuum_fit=false.')
    end

    chunk_size=get_effective_chunk_size(options,numel(list));

    [myflag,finaldata,fullwave,fullint,jd,signal_to_noise]=build_final_data_from_reduced_files(list,objectname,options,teff,logg,synthetic_vsini,sigm,wstart,wstop,shift,firstorder,lastorder,step,processing_directory,supplementary_directory,chunk_size);

    if myflag==0
        return
    end

else
    disp('previous final_data loaded')
    loaded_final_data=load(test1(1).name);
    finaldata=get_loaded_variable(loaded_final_data,'finaldata');
    if isfield(finaldata,'fullwave'), fullwave=finaldata.fullwave; end
    if isfield(finaldata,'fullint'), fullint=finaldata.fullint; end
    if isfield(finaldata,'jd'), jd=finaldata.jd; else, jd=NaN; end
    if isfield(finaldata,'signal_to_noise'), signal_to_noise=finaldata.signal_to_noise; end
end

%Guard line-profile requirements
if (options.calculate_lsd_profile==1 || options.radial_velocity_measurement==1 || options.FAMIAS_output==1) && options.order_merge~=1
    error('Regularised LSD profiles, radial velocity measurement, and FAMIAS output require options.order_merge=true.')
end

%OPTION:Regularised LSD profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.calculate_lsd_profile==1 || options.radial_velocity_measurement==1 || options.FAMIAS_output==1

    fprintf('Calculating regularised LSD profiles\n')

    if ~exist('fullwave','var') || ~exist('fullint','var')
        error('Regularised LSD profile calculation requires merged finaldata.fullwave and finaldata.fullint.')
    end

    %Find NaNs and set to continuum
    G=isnan(fullint);
    fullint(G)=1;

    %Create synthetic line mask for LSD
    lsd_line_mask=make_lsd_line_mask(teff,logg,synthetic_vsini,sigm,wstart,wstop,weak_limit1);

    %Set regularised LSD options
    regularised_lsd_options=get_regularised_lsd_options(options);

    %Calculate regularised LSD profiles
    regularised_lsd_profile=regularised_lsd_profile_from_spectra(fullwave,fullint,lsd_line_mask,jd,regularised_lsd_options);

    %Store product metadata
    regularised_lsd_profile.profile_type='LSD';
    regularised_lsd_profile.method='regularised_inverse_problem';
    regularised_lsd_profile.normalisation=regularised_lsd_options.normalisation;
    regularised_lsd_profile.objectname=objectname;
    regularised_lsd_profile.teff=teff;
    regularised_lsd_profile.logg=logg;
    regularised_lsd_profile.input_vsini=input_vsini;
    regularised_lsd_profile.synthetic_vsini=synthetic_vsini;

    finaldata.regularised_lsd_profile=regularised_lsd_profile;

end

%OPTION:Radial velocity and line broadening from regularised LSD profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.radial_velocity_measurement==1

    fprintf('Measuring radial velocities and vsini from regularised LSD profiles\n')

    if ~exist('regularised_lsd_profile','var')
        if isfield(finaldata,'regularised_lsd_profile')
            regularised_lsd_profile=finaldata.regularised_lsd_profile;
        else
            error('Missing regularised LSD profile in finaldata.')
        end
    end

    [radial_velocity,radial_velocity_error,systemic_velocity_direct_rv,systemic_velocity_gauss_rv,vsini,vsini_error,vbroad,vbroad_error,vsini_fit]=radial_velocity_calculation(regularised_lsd_profile);
    finaldata.radial_velocity=radial_velocity;
    finaldata.radial_velocity_error=radial_velocity_error;
    finaldata.systemic_velocity_direct=systemic_velocity_direct_rv;
    finaldata.systemic_velocity_gauss=systemic_velocity_gauss_rv;
    finaldata.vbroad=vbroad;
    finaldata.vbroad_error=vbroad_error;
    finaldata.vsini=vsini;
    finaldata.vsini_error=vsini_error;
    finaldata.vbroad_interpretation= ...
        ['Rotation-equivalent observed LSD width including instrumental ' ...
        'broadening and unresolved macroturbulence.'];
    finaldata.vsini_interpretation= ...
        ['Rotation-equivalent LSD width after subtracting the calibrated ' ...
        'instrumental variance; macroturbulence is not removed, so it is ' ...
        'not a pure projected rotational velocity.'];
    finaldata.vsini_instrumental_fwhm=vsini_fit.instrumental_fwhm;
    finaldata.vsini_instrumental_calibration_file= ...
        vsini_fit.instrumental_calibration_file;
    finaldata.rv_vsini_moment_diagnostics=vsini_fit;

    finaldata.single_vbroad=median(vbroad,'omitnan');
    finaldata.single_vsini=median(vsini,'omitnan');
    number_valid_vsini=sum(isfinite(vsini) & isfinite(vsini_error));
    if number_valid_vsini>0
        finaldata.single_vsini_error= ...
            median(vsini_error(isfinite(vsini) & isfinite(vsini_error)), ...
            'omitnan')/sqrt(number_valid_vsini);
    else
        finaldata.single_vsini_error=NaN;
    end
    number_valid_vbroad=sum(isfinite(vbroad) & isfinite(vbroad_error));
    if number_valid_vbroad>0
        finaldata.single_vbroad_error= ...
            median(vbroad_error(isfinite(vbroad) & isfinite(vbroad_error)), ...
            'omitnan')/sqrt(number_valid_vbroad);
    else
        finaldata.single_vbroad_error=NaN;
    end
    finaldata.single_vbroad_std=std(vbroad,'omitnan');
    finaldata.single_vsini_std=std(vsini,'omitnan');

    %Store LSD profiles in the stellar rest frame. The absolute systemic
    %velocity is retained separately, while the profile velocity axis and
    %coordinate-dependent diagnostics are shifted so the representative
    %profile is centred at zero.
    regularised_lsd_profile.velocity=regularised_lsd_profile.velocity- ...
        systemic_velocity_direct_rv;
    regularised_lsd_profile.systemic_velocity_removed= ...
        systemic_velocity_direct_rv;
    regularised_lsd_profile.systemic_velocity_reference='first moment of median LSD profile';
    regularised_lsd_profile.velocity_convention= ...
        'stellar rest-frame velocity; systemic velocity removed';
    finaldata.regularised_lsd_profile=regularised_lsd_profile;

    vsini_fit.velocity=vsini_fit.velocity-systemic_velocity_direct_rv;
    vsini_fit.line_window=vsini_fit.line_window-systemic_velocity_direct_rv;
    vsini_fit.profile_centre=vsini_fit.profile_centre-systemic_velocity_direct_rv;
    vsini_fit.profile_centre_centroid= ...
        vsini_fit.profile_centre_centroid-systemic_velocity_direct_rv;
    vsini_fit.moment_window=vsini_fit.moment_window-systemic_velocity_direct_rv;
    if isfield(vsini_fit,'representative_moments') && ...
            isfield(vsini_fit.representative_moments,'first_moment')
        vsini_fit.representative_moments.first_moment= ...
            vsini_fit.representative_moments.first_moment- ...
            systemic_velocity_direct_rv;
        if isfield(vsini_fit.representative_moments,'window')
            vsini_fit.representative_moments.window= ...
                vsini_fit.representative_moments.window- ...
                systemic_velocity_direct_rv;
        end
    end
    finaldata.rv_vsini_moment_diagnostics=vsini_fit;

end

%Plot the stored regularised LSD profiles. When RV measurement is enabled,
%the velocity axis above has already been placed in the stellar rest frame.
if exist('regularised_lsd_profile','var') && options.suppress_figures==0
    figure
    plot(regularised_lsd_profile.velocity,regularised_lsd_profile.intensity)
    title(['Regularised LSD profiles ' get_object_display_name(objectname)])
    xlabel('Velocity (km s^{-1})')
    ylabel('Regularised LSD intensity')
    grid on
    box on
end

%FINALISING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.order_merge==1
    if options.suppress_figures==0
        figure
        plot(fullwave,fullint)
        title(['Full Spectra ' get_object_display_name(objectname)])
    end
else
    if options.suppress_figures==0
        figure
        if isfield(finaldata,'wave87') && isfield(finaldata,'int87')
            plot(finaldata.wave87,finaldata.int87)
            title(['Sample Spectrum, order 87, ' get_object_display_name(objectname)])
        else
            warning(['Cannot plot sample order 87 for ' objectname ' because finaldata.wave87 or finaldata.int87 is missing.']);
        end
    end
end

%Plot radial velocities
if options.radial_velocity_measurement==1 && any(~isnan(radial_velocity(:)))
    if options.suppress_figures==0
        figure
        errorbar(jd,radial_velocity,make_error_vector(radial_velocity_error,jd),'xb')
        legend({'Regularised LSD'},'Location','best')
        title(['Radial Velocity ' get_object_display_name(objectname)])
        xlabel('JD')
        ylabel('Radial velocity')
        grid on
        box on
    end
end

%Plot measured line broadening
if options.radial_velocity_measurement==1
    if any(~isnan(vsini(:)))
        if options.suppress_figures==0
            figure
            errorbar(jd,vsini,make_error_vector(vsini_error,jd),'x')
            legend({'Regularised LSD'},'Location','best')
            title(['vsini ' get_object_display_name(objectname)])
            xlabel('JD')
            ylabel('vsini (km s^{-1})')
            grid on
            box on
        end
    end
end

%Mark processing complete
myflag=1;

%Store signal-to-noise values
if isstruct(finaldata)
    finaldata.signal_to_noise=signal_to_noise;
end

%Ensure output exists
if exist('finaldata','var')==0
    finaldata=[];
end

%Keep diagnostic and fit internals out of the main final-data product.
[finaldata,new_diagnostics]=separate_finaldata_diagnostics(finaldata);
diagnostics_file=fullfile(supplementary_directory, ...
    ['diagnostics_' objectname '.mat']);
diagnostics=merge_saved_diagnostics(diagnostics_file,new_diagnostics);
save(diagnostics_file,'diagnostics','-v7.3')

fprintf('DATA PROCESSED\n')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [myflag,finaldata,fullwave,fullint,jd,signal_to_noise]=build_final_data_from_reduced_files(list,objectname,options,teff,logg,synthetic_vsini,sigm,wstart,wstop,shift,firstorder,lastorder,step,processing_directory,supplementary_directory,chunk_size)

myflag=1;
finaldata=[];
fullwave=[];
fullint=[];
jd=[];
signal_to_noise=[];

number_observations=numel(list);
frac_thresh=0.20;

fprintf('Pre-cleaning reduced spectra before post-reduction\n')

[all_minwaves,all_maxwaves]=measure_order_wavelength_limits(list,firstorder,lastorder);
[order_edge_variation,order_edge_report]=find_order_edge_rejections( ...
    all_minwaves,all_maxwaves,frac_thresh);

if any(order_edge_variation)
    save(fullfile(supplementary_directory,'post_reduction_order_edge_coverage.mat'), ...
        'order_edge_report','order_edge_variation')
    fprintf(['Preserving %.0f spectra with differing order-edge coverage; ' ...
        'uncovered pixels will have zero weight.\n'],sum(order_edge_variation))
else
    fprintf('All spectra have consistent order-edge coverage\n')
end

minwave=min(all_minwaves,[],1,'omitnan');
maxwave=max(all_maxwaves,[],1,'omitnan');

data_grid=define_union_order_grids(minwave,maxwave,firstorder,lastorder,step);

[datacf,continuum_fit_source]=load_continuum_fit_for_processing(options,objectname);

order_medians=compute_global_order_medians_sampled(list,objectname,options,firstorder,lastorder,step,data_grid,500);

chunk_directory=fullfile(supplementary_directory,'chunked_post_reduction_files');
if ~isfolder(chunk_directory), mkdir(chunk_directory); end
delete_existing_chunk_files(chunk_directory)

number_chunks=ceil(number_observations/chunk_size);
if number_chunks>1
    fprintf(['Processing %.0f reduced spectra in %.0f blocks of up to ' ...
        '%.0f spectra\n'],number_observations,number_chunks,chunk_size)
else
    fprintf('Processing %.0f reduced spectra\n',number_observations)
end

row_start=1;
cosmic_ray_cleanup_chunks=cell(number_chunks,1);

for chunk_index=1:number_chunks

    chunk_start=(chunk_index-1)*chunk_size+1;
    chunk_end=min(chunk_index*chunk_size,number_observations);
    chunk_list=list(chunk_start:chunk_end);

    if number_chunks>1
        fprintf('Processing block %.0f of %.0f: spectra %.0f to %.0f\n', ...
            chunk_index,number_chunks,chunk_start,chunk_end)
    end

    [chunk_fullwave,chunk_fullint,chunk_jd,chunk_signal_to_noise,datacf,chunk_cosmic_ray_cleanup]=process_reduced_file_chunk(chunk_list,objectname,options,teff,logg,synthetic_vsini,sigm,wstart,wstop,shift,firstorder,lastorder,step,data_grid,order_medians,datacf,continuum_fit_source,processing_directory);
    cosmic_ray_cleanup_chunks{chunk_index}=chunk_cosmic_ray_cleanup;
    continuum_fit_source='loaded';

    if isempty(fullwave)
        fullwave=chunk_fullwave;
        fullint=nan(number_observations,size(chunk_fullint,2));
        jd=nan(number_observations,1);
        signal_to_noise=nan(number_observations,1);
    else
        if numel(chunk_fullwave)~=numel(fullwave) || any(abs(chunk_fullwave(:)-fullwave(:))>1e-10)
            error('Merged wavelength grid changed between chunks. Check combining_orders_2Dint and common order grids.')
        end
    end

    row_end=row_start+size(chunk_fullint,1)-1;
    fullint(row_start:row_end,:)=chunk_fullint;
    jd(row_start:row_end)=chunk_jd;
    signal_to_noise(row_start:row_end)=chunk_signal_to_noise;

    chunkdata=struct();
    chunkdata.fullwave=chunk_fullwave;
    chunkdata.fullint=chunk_fullint;
    chunkdata.jd=chunk_jd;
    chunkdata.signal_to_noise=chunk_signal_to_noise;
    chunkdata.source_files={chunk_list.name};
    chunkdata.cosmic_ray_cleanup=chunk_cosmic_ray_cleanup;
    save(fullfile(chunk_directory,['chunk_' sprintf('%05d',chunk_index) '.mat']),'chunkdata','-v7.3')

    row_start=row_end+1;

end

save(fullfile(supplementary_directory,'post_reduction_cosmic_ray_cleanup.mat'), ...
    'cosmic_ray_cleanup_chunks','-v7.3')

[fullint,merged_pixel_cleanup.before_suppnet]= ...
    mask_nonphysical_merged_pixels(fullwave,fullint);
[fullint,merged_pixel_cleanup.positive_spikes_before_suppnet]= ...
    remove_isolated_merged_cosmic_ray_spikes(fullwave,fullint);

merged_reject=find_merged_spectrum_rejections(fullwave,fullint);

if any(merged_reject)

    move_rejected_reduced_files(list,find(merged_reject), ...
        fullfile(fileparts(processing_directory),'removed_reduced_files'), ...
        fullfile(processing_directory,'removed_reduced_files'));
    save_cleanup_report(fullfile(supplementary_directory, ...
        'post_reduction_merged_spectrum_cleanup.mat'),list,struct(),merged_reject)

    fprintf('A total of %.0f observations have been removed from the merged spectra\n',sum(merged_reject))

    keep_merged=~merged_reject;

    if ~any(keep_merged)
        error('All merged spectra were rejected. Check post-reduction cleanup thresholds.')
    end

    fullint=fullint(keep_merged,:);
    jd=jd(keep_merged);
    signal_to_noise=signal_to_noise(keep_merged);

else

    fprintf('A total of 0 observations have been removed from the merged spectra\n')

end

if options.automatic_continuum_fit==1
    fprintf('Applying SUPPNet to the full merged spectrum\n')
    [fullint,merged_suppnet]=apply_suppnet_merged_continuum( ...
        fullwave,fullint,objectname);
    [fullint,merged_pixel_cleanup.after_suppnet]= ...
        mask_nonphysical_merged_pixels(fullwave,fullint);
    finaldata.merged_suppnet=merged_suppnet;
    if isfield(datacf,'automatic_fit_method')
        finaldata.automatic_fit_method=datacf.automatic_fit_method;
    end
    if isfield(datacf,'automatic_fit_summary')
        finaldata.automatic_fit_summary=datacf.automatic_fit_summary;
    end
    if isfield(datacf,'automatic_order_scales')
        finaldata.automatic_order_scales=datacf.automatic_order_scales;
    end
elseif options.template_continuum_fit==1
    if isfield(datacf,'template_fit_method')
        finaldata.template_fit_method=datacf.template_fit_method;
    end
    if isfield(datacf,'template_fit_summary')
        finaldata.template_fit_summary=datacf.template_fit_summary;
    end
end
save(fullfile(supplementary_directory, ...
    'post_reduction_merged_pixel_cleanup.mat'),'merged_pixel_cleanup')

finaldata.fullwave=fullwave';
finaldata.fullint=fullint;
finaldata.jd=jd;
finaldata.extended_wavelength_range=options.extended_wavelength_range;
if options.extended_wavelength_range
    finaldata.requested_wavelength_limits=[3750 8950];
else
    finaldata.requested_wavelength_limits=[3869 6866];
end

end

function [all_minwaves,all_maxwaves]=measure_order_wavelength_limits(list,firstorder,lastorder)

number_observations=numel(list);
number_orders=lastorder-firstorder+1;
all_minwaves=nan(number_observations,number_orders);
all_maxwaves=nan(number_observations,number_orders);

for obsnum=1:number_observations

    reduced_frame=load_primary_variable_from_mat_file(get_list_file_path(list(obsnum)));
    cntt=0;

    for ord=firstorder:lastorder

        cntt=cntt+1;
        order_field_name=['order_' num2str(ord)];

        if ~isfield(reduced_frame,order_field_name)
            continue
        end

        wave=reduced_frame.(order_field_name)(:,1);
        all_minwaves(obsnum,cntt)=min(wave);
        all_maxwaves(obsnum,cntt)=max(wave);

    end

    clear reduced_frame wave

end

end

function [reject_file,cleanup_report]=find_order_edge_rejections(all_minwaves,all_maxwaves,frac_thresh)

number_observations=size(all_minwaves,1);
number_orders=size(all_minwaves,2);
bad_min=false(number_observations,number_orders);
bad_max=false(number_observations,number_orders);
threshold_per_order=nan(1,number_orders);
median_min=nan(1,number_orders);
median_max=nan(1,number_orders);

for order_index=1:number_orders

    this_min=all_minwaves(:,order_index);
    this_max=all_maxwaves(:,order_index);
    med_min=median(this_min,'omitnan');
    med_max=median(this_max,'omitnan');
    span_med=med_max-med_min;
    thresh=frac_thresh*span_med;

    threshold_per_order(order_index)=thresh;
    median_min(order_index)=med_min;
    median_max(order_index)=med_max;

    bad_min(:,order_index)=(this_min-med_min)>thresh;
    bad_max(:,order_index)=(med_max-this_max)>thresh;

end

bad_order=bad_min | bad_max;
reject_file=any(bad_order,2);

cleanup_report=struct();
cleanup_report.bad_min=bad_min;
cleanup_report.bad_max=bad_max;
cleanup_report.threshold_per_order=threshold_per_order;
cleanup_report.median_min=median_min;
cleanup_report.median_max=median_max;

end

function data_grid=define_union_order_grids(minwave,maxwave,firstorder,lastorder,step)

data_grid=struct();

for ord=firstorder:lastorder

    order_index=ord-firstorder+1;
    wave_field_name=['wave' num2str(ord)];

    if isnan(minwave(order_index)) || isnan(maxwave(order_index)) || minwave(order_index)+2>=maxwave(order_index)-2
        data_grid.(wave_field_name)=[];
    else
        data_grid.(wave_field_name)=rounding(minwave(order_index)+2,step):step:rounding(maxwave(order_index)-2,step);
    end

end

end

function [datacf,continuum_fit_source]=load_continuum_fit_for_processing(options,objectname)

datacf=[];
continuum_fit_source='';
customdatacf=fullfile(pwd,'supplementary_data',['datacf_' objectname]);

fprintf('CONTINUUM FITTING\n')

if options.manual_continuum_fit==1
    continuum_fit_source='manual';
elseif options.automatic_continuum_fit==1
    continuum_fit_source='automatic';
elseif options.template_continuum_fit==1
    continuum_fit_source='template';
elseif isfile([customdatacf '.mat']) || isfile(customdatacf)
    loaded_continuum_fit=load(customdatacf);
    datacf=get_loaded_variable(loaded_continuum_fit,'datacf');
    continuum_fit_source='loaded';
    fprintf('Loaded previous manual continuum fit\n')
else
    error('No previous manual continuum fit found, and no alternative continuum option selected.')
end

end

function order_medians=compute_global_order_medians_sampled(list,objectname,options,firstorder,lastorder,step,data_grid,maximum_reference_spectra)

order_medians=struct();

number_files=numel(list);
number_reference_spectra=min(maximum_reference_spectra,number_files);
reference_indices=unique(round(linspace(1,number_files,number_reference_spectra)));
reference_list=list(reference_indices);
number_reference_spectra=numel(reference_list);

order_matrices=struct();

for ord=firstorder:lastorder

    wave_field_name=['wave' num2str(ord)];
    matrix_field_name=['order' num2str(ord)];

    if ~isfield(data_grid,wave_field_name) || isempty(data_grid.(wave_field_name))
        order_matrices.(matrix_field_name)=[];
        continue
    end

    wave=data_grid.(wave_field_name);

    if numel(wave)<=1000
        order_matrices.(matrix_field_name)=[];
        continue
    end

    wave_for_order=wave(50:end-50);
    order_matrices.(matrix_field_name)=nan(number_reference_spectra,numel(wave_for_order));

end

for obsnum=1:number_reference_spectra

    reduced_file_name=get_list_file_path(reference_list(obsnum));
    file=get_reduced_file_base_name(reduced_file_name);
    reduced_frame=load_primary_variable_from_mat_file(reduced_file_name);

    for ord=firstorder:lastorder

        order_field_name=['order_' num2str(ord)];
        wave_field_name=['wave' num2str(ord)];
        matrix_field_name=['order' num2str(ord)];

        if ~isfield(reduced_frame,order_field_name)
            continue
        end

        if ~isfield(data_grid,wave_field_name) || isempty(data_grid.(wave_field_name))
            continue
        end

        if ~isfield(order_matrices,matrix_field_name) || isempty(order_matrices.(matrix_field_name))
            continue
        end

        wave=data_grid.(wave_field_name);
        wave_for_order=wave(50:end-50);

        if options.apply_barycentric_correction==1
            [tmpwave,~,~]=get_barycentric_corrected_wave(reduced_frame,objectname,file,ord);
        else
            tmpwave=reduced_frame.(order_field_name)(:,1);
        end

        int=reduced_frame.(order_field_name)(:,2);
        order_matrices.(matrix_field_name)(obsnum,:)= ...
            interpolate_without_extrapolation(tmpwave,int,wave_for_order,NaN);

    end

    clear reduced_frame

end

for ord=firstorder:lastorder

    wave_field_name=['wave' num2str(ord)];
    median_field_name=['median' num2str(ord)];
    matrix_field_name=['order' num2str(ord)];

    if ~isfield(data_grid,wave_field_name) || isempty(data_grid.(wave_field_name))
        order_medians.(median_field_name)=[];
        continue
    end

    wave=data_grid.(wave_field_name);

    if numel(wave)<=1000
        order_medians.(median_field_name)=ones(1,max(numel(wave)-140,0));
        continue
    end

    order_matrix=order_matrices.(matrix_field_name);

    order_matrix(order_matrix>2)=2;
    intx=medfilt1(order_matrix',71)';
    intxx=intx(:,71:end-70);
    medianspec=median(intxx,1,'omitnan');

    bad_median=~isfinite(medianspec) | medianspec==0;

    if any(bad_median)
        medianspec(bad_median)=1;
    end

    order_medians.(median_field_name)=medianspec;

    clear order_matrix intx intxx medianspec

end

end

function [fullwave,fullint,jd,signal_to_noise,datacf,cosmic_ray_cleanup]=process_reduced_file_chunk(list,objectname,options,teff,logg,synthetic_vsini,sigm,wstart,wstop,shift,firstorder,lastorder,step,data_grid,order_medians,datacf,continuum_fit_source,processing_directory)

data=struct();
jd=nan(numel(list),1);
signal_to_noise=nan(numel(list),1);
allbcorr=nan(numel(list),1); %#ok<NASGU>

%Preallocate the union grids. Individual spectra retain NaN intensity and
%zero weight wherever their extracted order does not provide data.
for ord=firstorder:lastorder
    wave_field_name=['wave' num2str(ord)];
    intensity_field_name=['int' num2str(ord)];
    weights_field_name=['weights' num2str(ord)];
    if ~isfield(data_grid,wave_field_name) || ...
            isempty(data_grid.(wave_field_name)), continue; end
    data.(wave_field_name)=data_grid.(wave_field_name);
    data.(intensity_field_name)=nan(numel(list), ...
        numel(data_grid.(wave_field_name)));
    data.(weights_field_name)=zeros(numel(list), ...
        numel(data_grid.(wave_field_name)));
end

for obsnum=1:numel(list)

    reduced_file_name=get_list_file_path(list(obsnum));
    file=get_reduced_file_base_name(reduced_file_name);
    reduced_frame=load_primary_variable_from_mat_file(reduced_file_name);
    processing_frame=load_processing_frame(processing_directory,file);
    bcorr=NaN;

    for ord=firstorder:lastorder

        order_field_name=['order_' num2str(ord)];
        wave_field_name=['wave' num2str(ord)];
        intensity_field_name=['int' num2str(ord)];
        weights_field_name=['weights' num2str(ord)];

        if ~isfield(reduced_frame,order_field_name), continue; end
        if ~isfield(data_grid,wave_field_name) || isempty(data_grid.(wave_field_name)), continue; end

        if options.apply_barycentric_correction==1
            [tmpwave,bcorr,reduced_frame]=get_barycentric_corrected_wave(reduced_frame,objectname,file,ord);
        else
            tmpwave=reduced_frame.(order_field_name)(:,1);
        end

        int=reduced_frame.(order_field_name)(:,2);
        cutpixel_order_field_name=['order_' num2str(ord)];
        normalised_flat_field_name=['order_nff_' num2str(ord)];
        flat_quality_field_name=['order_flat_quality_' num2str(ord)];
        smoothed_flat_field_name=['order_flat_smoothed_' num2str(ord)];
        extraction_quality_field_name=['order_extraction_quality_' num2str(ord)];

        if ~isfield(processing_frame,'order_cutpixels') || ~isfield(processing_frame.order_cutpixels,cutpixel_order_field_name), error(['Missing order_cutpixels.' cutpixel_order_field_name ' in processing file for ' file '.']); end
        if ~isfield(processing_frame,normalised_flat_field_name), error(['Missing ' normalised_flat_field_name ' in processing file for ' file '.']); end

        iind1=processing_frame.order_cutpixels.(cutpixel_order_field_name).ind1;
        iind2=processing_frame.order_cutpixels.(cutpixel_order_field_name).ind2;
        signal_weight=processing_frame.(normalised_flat_field_name)(iind1:iind2);
        if isfield(processing_frame,flat_quality_field_name)
            flat_quality=processing_frame.(flat_quality_field_name)(iind1:iind2);
        elseif isfield(processing_frame,smoothed_flat_field_name)
            flat_quality=calculate_flat_quality( ...
                processing_frame.(smoothed_flat_field_name)(iind1:iind2));
        else
            error(['Missing flat-quality information for ' file ...
                ', order ' num2str(ord) '.'])
        end
        if isfield(processing_frame,extraction_quality_field_name)
            extraction_quality=processing_frame.( ...
                extraction_quality_field_name)(iind1:iind2);
        else
            extraction_quality=calculate_extraction_quality_from_processing( ...
                processing_frame,ord,iind1,iind2,int);
        end
        weights=max(signal_weight,0).*flat_quality.^2.* ...
            extraction_quality.^2;
        int(flat_quality<=0 | extraction_quality<=0)=NaN;

        if numel(tmpwave)~=numel(int) || numel(tmpwave)~=numel(weights)
            error(['Wavelength, intensity and flat-weight lengths differ for ' ...
                file ', order ' num2str(ord) '.'])
        end

        interpolated_intensity=interpolate_without_extrapolation( ...
            tmpwave,int,data_grid.(wave_field_name),NaN);
        interpolated_weights=interpolate_without_extrapolation( ...
            tmpwave,weights,data_grid.(wave_field_name),0);
        invalid_pixels=~isfinite(interpolated_intensity) | ...
            ~isfinite(interpolated_weights) | interpolated_weights<=0;
        interpolated_intensity(invalid_pixels)=NaN;
        interpolated_weights(invalid_pixels)=0;
        data.(intensity_field_name)(obsnum,:)=interpolated_intensity;
        data.(weights_field_name)(obsnum,:)=interpolated_weights;

    end

    if isfield(reduced_frame,'jd'), jd(obsnum)=reduced_frame.jd; end
    if isfield(reduced_frame,'signal_to_noise'), signal_to_noise(obsnum)=reduced_frame.signal_to_noise; end
    allbcorr(obsnum)=bcorr; %#ok<NASGU>

    if isfield(reduced_frame,'bcorr_was_updated') && reduced_frame.bcorr_was_updated
        reduced_frame=rmfield(reduced_frame,'bcorr_was_updated');
        save_reduced_frame_with_original_variable_name(reduced_file_name,file,reduced_frame);
    end

    clear reduced_frame processing_frame

end

fprintf('Normalising all orders\n')

data2=normalise_orders_with_global_references(data,order_medians,firstorder,lastorder);
clear data
data2=apply_extended_order_limits(data2,firstorder,lastorder, ...
    options.extended_wavelength_range);

if strcmp(continuum_fit_source,'manual')
    fprintf('Create new manual continuum fit.\n')
    fprintf('Opening program for manual fit\n')
    fprintf('Try to keep the function smooth and fit the top sides of the line profiles.\n')
    datacf=manual_continuum_fitting(teff,logg,synthetic_vsini,sigm,wstart,wstop,shift,firstorder,lastorder,data2,objectname);
elseif strcmp(continuum_fit_source,'automatic')
    fprintf('Creating synthetic-guided coarse continuum fit\n')
    datacf=prepare_coarse_continuum_fit(data2,objectname,teff,logg, ...
        synthetic_vsini,sigm,wstart,wstop,shift,firstorder,lastorder);
elseif strcmp(continuum_fit_source,'template')
    fprintf('Interpolating the SpT continuum-template grid\n')
    datacf=prepare_spt_template_continuum_fit(data2,objectname,teff, ...
        firstorder,lastorder,options.extended_wavelength_range);
end

data2=apply_continuum_fit_to_orders(data2,datacf,firstorder,lastorder);
[data2,cosmic_ray_cleanup]=remove_isolated_cosmic_ray_spikes( ...
    data2,firstorder,lastorder);
if options.automatic_continuum_fit==1
    [data2,datacf]=harmonise_continuum_orders(data2,datacf, ...
        firstorder,lastorder,step);
end
data2=mask_nonphysical_order_pixels(data2,firstorder,lastorder);
[data2,effective_firstorder,effective_lastorder]=prune_empty_orders( ...
    data2,firstorder,lastorder);

%Run the optional pre-merge order RV diagnostic without changing data2.
if isfield(options,'order_rv_diagnostic') && options.order_rv_diagnostic
    try
        evaluate_order_by_order_rv(data2,jd,objectname,teff,logg, ...
            synthetic_vsini,sigm,wstart,wstop);
    catch order_rv_diagnostic_error
        warning(['Optional order-by-order RV diagnostic was skipped: ' ...
            order_rv_diagnostic_error.message])
    end
end

fprintf('Merging orders\n')
[fullwave,fullint]=combining_orders_2Dint(data2,step, ...
    effective_firstorder,effective_lastorder);

end

function flat_quality=calculate_flat_quality(smoothed_flat)

%Convert normalised flat response into an inverse-variance edge taper
flat_quality_floor=0.2;
flat_quality_full=0.5;
flat_quality=zeros(size(smoothed_flat));
valid_flat=isfinite(smoothed_flat) & smoothed_flat>flat_quality_floor;
flat_quality(valid_flat)=(smoothed_flat(valid_flat)-flat_quality_floor)./ ...
    (flat_quality_full-flat_quality_floor);
flat_quality=min(max(flat_quality,0),1);

end

function extraction_quality=calculate_extraction_quality_from_processing( ...
    processing_frame,order_number,first_pixel,last_pixel,reduced_intensity)

%Reconstruct extraction uncertainty for processing files made before flags
source_field=['order_nff_' num2str(order_number)];
background_field=['order_bkgd_smoothed_' num2str(order_number)];
flat_field=['order_flat_smoothed_' num2str(order_number)];

source_counts=processing_frame.(source_field)(first_pixel:last_pixel);
background_counts=processing_frame.(background_field)(first_pixel:last_pixel);
smoothed_flat=processing_frame.(flat_field)(first_pixel:last_pixel);
xmed_values=(source_counts-background_counts)./ ...
    (smoothed_flat.*reduced_intensity);
valid_xmed=isfinite(xmed_values) & xmed_values>0 & ...
    isfinite(smoothed_flat) & smoothed_flat>=0.5 & ...
    isfinite(reduced_intensity) & abs(reduced_intensity)>0.05;

if nnz(valid_xmed)<20
    extraction_quality=ones(size(reduced_intensity));
    return
end

xmed=median(xmed_values(valid_xmed),'omitnan');
normalised_uncertainty=sqrt(max(source_counts,0)+ ...
    abs(background_counts))./(abs(xmed).*smoothed_flat);
full_quality_limit=0.10;
zero_quality_limit=0.20;
extraction_quality=(zero_quality_limit-normalised_uncertainty)./ ...
    (zero_quality_limit-full_quality_limit);
extraction_quality=min(max(extraction_quality,0),1);
extraction_quality(~isfinite(normalised_uncertainty))=0;

end

function data2=normalise_orders_with_global_references(data,order_medians,firstorder,lastorder)

data2=struct();

for ord=firstorder:lastorder

    small_flag=0;
    wave_field_name=['wave' num2str(ord)];
    intensity_field_name=['int' num2str(ord)];
    weights_field_name=['weights' num2str(ord)];
    median_field_name=['median' num2str(ord)];

    if ~isfield(data,wave_field_name) || ~isfield(data,intensity_field_name), continue; end

    int=data.(intensity_field_name);
    wave=data.(wave_field_name);

    if size(int,2)>1000
        int=int(:,50:end-50);
        wave=wave(50:end-50);
        order_weights=data.(weights_field_name)(:,50:end-50);
        small_flag=1;
    else
        order_weights=data.(weights_field_name);
    end

    int(int>2)=2;
    intx=medfilt1(int',71)';
    intxx=intx(:,71:end-70);
    wavex=wave(71:end-70);
    fit_weights=order_weights(:,71:end-70);

    if isfield(order_medians,median_field_name) && numel(order_medians.(median_field_name))==numel(wavex)
        medianspec=order_medians.(median_field_name);
    else
        medianspec=median(intxx,1,'omitnan');
    end

    bad_median=~isfinite(medianspec) | medianspec==0;
    if any(bad_median)
        medianspec(bad_median)=1;
    end

    if small_flag==1
        output_weights=data.(weights_field_name)(:,50:end-50);
    else
        output_weights=data.(weights_field_name);
    end
    normalised_intensity=nan(size(int));
    for ind1=1:size(int,1)

        temp=intxx(ind1,:);
        temp_m=temp./medianspec;
        telluric_contamination=wavex>=7590 & wavex<=7700;
        valid_fit=isfinite(wavex) & isfinite(temp_m) & temp_m>0 & ...
            fit_weights(ind1,:)>0 & ~telluric_contamination;
        if nnz(valid_fit)>=6
            ratio_limits=prctile(temp_m(valid_fit),[10 90]);
            robust_fit=valid_fit & temp_m>=ratio_limits(1) & ...
                temp_m<=ratio_limits(2);
            if nnz(robust_fit)<6, robust_fit=valid_fit; end
            warning off
            fit_values=temp_m(robust_fit);
            [p,S,mu]=polyfit(wavex(robust_fit),log(fit_values),3);
            warning on
            thefit=exp(polyval(p,wave,S,mu));
            fit_wave_min=min(wavex(robust_fit));
            fit_wave_max=max(wavex(robust_fit));
            thefit(wave<fit_wave_min)=exp(polyval(p,fit_wave_min,S,mu));
            thefit(wave>fit_wave_max)=exp(polyval(p,fit_wave_max,S,mu));

            %Keep the smooth correction within the range supported by the
            %observed order ratio. This does not clip the spectrum itself.
            fit_limits=prctile(fit_values,[5 95]);
            lower_fit_limit=max(0.5*fit_limits(1),eps);
            upper_fit_limit=2*fit_limits(2);
            thefit=min(max(thefit,lower_fit_limit),upper_fit_limit);
        else
            thefit=ones(size(wave));
        end
        bad_fit=~isfinite(thefit) | thefit<=0;
        thefit(bad_fit)=NaN;
        into=int(ind1,:)./thefit;
        output_weights(ind1,bad_fit | ~isfinite(into))=0;
        normalised_intensity(ind1,:)=into;

    end

    data2.(intensity_field_name)=normalised_intensity;

    data2.(wave_field_name)=wave;

    data2.(weights_field_name)=output_weights;

end

end

function interpolated_values=interpolate_without_extrapolation( ...
    wavelength,values,output_wavelength,outside_value)

wavelength=wavelength(:);
values=values(:);
valid_wavelength=isfinite(wavelength);
wavelength=wavelength(valid_wavelength);
values=values(valid_wavelength);

if numel(wavelength)<2
    interpolated_values=outside_value*ones(size(output_wavelength));
    return
end

[wavelength,sort_index]=sort(wavelength);
values=values(sort_index);
[wavelength,unique_index]=unique(wavelength,'stable');
values=values(unique_index);

interpolated_values=outside_value*ones(size(output_wavelength));
valid_values=isfinite(values);
segment_start=find(valid_values & [true;~valid_values(1:end-1)]);
segment_stop=find(valid_values & [~valid_values(2:end);true]);

for segment_index=1:numel(segment_start)
    first_index=segment_start(segment_index);
    last_index=segment_stop(segment_index);
    if last_index-first_index+1<2, continue; end
    segment_wave=wavelength(first_index:last_index);
    segment_values=values(first_index:last_index);
    inside=output_wavelength>=segment_wave(1) & ...
        output_wavelength<=segment_wave(end);
    interpolated_values(inside)=interp1(segment_wave,segment_values, ...
        output_wavelength(inside),'linear');
end

end

function data2=apply_continuum_fit_to_orders(data2,datacf,firstorder,lastorder)

warning off
for ord=firstorder:lastorder

    wave_field_name=['wave' num2str(ord)];
    intensity_field_name=['int' num2str(ord)];
    weights_field_name=['weights' num2str(ord)];
    continuum_intensity_field_name=['int' num2str(ord)];
    continuum_wave_field_name=['w' num2str(ord)];

    if ~isfield(data2,wave_field_name) || ~isfield(data2,intensity_field_name), continue; end

    if ~isfield(datacf,continuum_intensity_field_name) || ~isfield(datacf,continuum_wave_field_name)
        warning(['Continuum fit is missing order ' num2str(ord) '. Using unity continuum for this order.'])
        data2.(intensity_field_name)=apply_continuum_vector_to_intensity_matrix(data2.(intensity_field_name),ones(1,size(data2.(intensity_field_name),2)),intensity_field_name);
        continue
    end

    wave=data2.(wave_field_name);
    fiti=datacf.(continuum_intensity_field_name);
    fitw=datacf.(continuum_wave_field_name);

    if isempty(fiti) || isempty(fitw)
        warning(['Continuum fit for order ' num2str(ord) ' is empty. Using unity continuum for this order.'])
        data2.(intensity_field_name)=apply_continuum_vector_to_intensity_matrix(data2.(intensity_field_name),ones(1,size(data2.(intensity_field_name),2)),intensity_field_name);
        continue
    end

    if wave(1)~=fitw(1) || numel(wave)~=numel(fitw)
        continuum_vector=interpolate_continuum_with_endpoint_hold( ...
            fitw,fiti,wave);
    else
        continuum_vector=fiti;
    end
    continuum_vector=continuum_vector(:)';

    invalid_continuum=~isfinite(continuum_vector) | continuum_vector<=0;
    continuum_vector(invalid_continuum)=NaN;
    data2.(intensity_field_name)=apply_continuum_vector_to_intensity_matrix( ...
        data2.(intensity_field_name),continuum_vector,intensity_field_name);
    if isfield(data2,weights_field_name)
        invalid_intensity=~isfinite(data2.(intensity_field_name));
        invalid_pixels=repmat(invalid_continuum, ...
            size(data2.(intensity_field_name),1),1) | invalid_intensity;
        data2.(weights_field_name)(invalid_pixels)=0;
    end

end
warning on

end

function continuum=interpolate_continuum_with_endpoint_hold( ...
    fit_wavelength,fit_intensity,output_wavelength)

fit_wavelength=fit_wavelength(:);
fit_intensity=fit_intensity(:);
valid_fit=isfinite(fit_wavelength) & isfinite(fit_intensity) & ...
    fit_intensity>0;
fit_wavelength=fit_wavelength(valid_fit);
fit_intensity=fit_intensity(valid_fit);

if numel(fit_wavelength)<2
    continuum=nan(size(output_wavelength));
    return
end

[fit_wavelength,sort_index]=sort(fit_wavelength);
fit_intensity=fit_intensity(sort_index);
[fit_wavelength,unique_index]=unique(fit_wavelength,'stable');
fit_intensity=fit_intensity(unique_index);

continuum=interp1(fit_wavelength,fit_intensity,output_wavelength, ...
    'linear',NaN);
continuum(output_wavelength<fit_wavelength(1))=fit_intensity(1);
continuum(output_wavelength>fit_wavelength(end))=fit_intensity(end);

end

function data2=apply_extended_order_limits( ...
    data2,firstorder,lastorder,extended_wavelength_range)
%Apply standard limits unless the manually reviewed extensions are enabled.

if extended_wavelength_range
    blue_wavelength_limit=3750;
    red_wavelength_limit=8950;
else
    blue_wavelength_limit=3869;
    red_wavelength_limit=6866;
end

for order_number=firstorder:lastorder
    wave_field_name=['wave' num2str(order_number)];
    intensity_field_name=['int' num2str(order_number)];
    weights_field_name=['weights' num2str(order_number)];
    if ~isfield(data2,wave_field_name) || ...
            ~isfield(data2,intensity_field_name) || ...
            ~isfield(data2,weights_field_name)
        continue
    end

    wave=data2.(wave_field_name);
    keep=wave>=blue_wavelength_limit & wave<=red_wavelength_limit;

    if nnz(keep)<2
        data2=rmfield(data2,{wave_field_name,intensity_field_name, ...
            weights_field_name});
        continue
    end

    data2.(wave_field_name)=wave(keep);
    data2.(intensity_field_name)=data2.(intensity_field_name)(:,keep);
    data2.(weights_field_name)=data2.(weights_field_name)(:,keep);
end

end

function [data2,effective_firstorder,effective_lastorder]= ...
    prune_empty_orders(data2,firstorder,lastorder)

effective_firstorder=firstorder;
effective_lastorder=lastorder;

for ind3=lastorder:-1:firstorder

    wave_field_name=['wave' num2str(ind3)];
    intensity_field_name=['int' num2str(ind3)];
    weights_field_name=['weights' num2str(ind3)];

    if ~isfield(data2,wave_field_name) || isempty(data2.(wave_field_name))
        if isfield(data2,wave_field_name), data2=rmfield(data2,wave_field_name); end
        if isfield(data2,intensity_field_name), data2=rmfield(data2,intensity_field_name); end
        if isfield(data2,weights_field_name), data2=rmfield(data2,weights_field_name); end
        if ind3==effective_lastorder
            effective_lastorder=ind3-1;
        end
    end

end

for order_number=firstorder:lastorder
    wave_field_name=['wave' num2str(order_number)];
    if isfield(data2,wave_field_name) && ~isempty(data2.(wave_field_name))
        effective_firstorder=order_number;
        break
    end
end

if effective_lastorder<effective_firstorder
    error('No valid orders remain after pruning empty order fields.')
end

end

function data2=mask_nonphysical_order_pixels(data2,firstorder,lastorder)

%Give non-physical negative order pixels zero weight before merging
minimum_normalised_intensity=-0.1;
number_masked=0;

for order_number=firstorder:lastorder
    intensity_field_name=['int' num2str(order_number)];
    weights_field_name=['weights' num2str(order_number)];
    if ~isfield(data2,intensity_field_name) || ...
            ~isfield(data2,weights_field_name)
        continue
    end
    intensity=data2.(intensity_field_name);
    bad_pixel=isfinite(intensity) & ...
        intensity<minimum_normalised_intensity;
    number_masked=number_masked+sum(bad_pixel(:));
    intensity(bad_pixel)=NaN;
    data2.(intensity_field_name)=intensity;
    data2.(weights_field_name)(bad_pixel)=0;
end

fprintf(['Masked %.0f non-physical order pixels locally; ' ...
    'the remaining wavelengths in each observation were retained\n'], ...
    number_masked)

end

function [fullint,cleanup]=mask_nonphysical_merged_pixels(fullwave,fullint)

%Retain each observation while removing non-physical negative pixels
cleanup.lower_limit=-0.1;

bad_pixel=isfinite(fullint) & fullint<cleanup.lower_limit;

cleanup.number_pixels=sum(bad_pixel(:));
cleanup.affected_spectra=find(any(bad_pixel,2));
cleanup.affected_wavelengths=find(any(bad_pixel,1));

if isempty(cleanup.affected_wavelengths)
    cleanup.minimum_wavelength=NaN;
    cleanup.maximum_wavelength=NaN;
else
    cleanup.minimum_wavelength=fullwave(cleanup.affected_wavelengths(1));
    cleanup.maximum_wavelength=fullwave(cleanup.affected_wavelengths(end));
end

fullint(bad_pixel)=NaN;

if cleanup.number_pixels>0
    fprintf(['Masked %.0f non-physical merged pixels in %.0f spectra ' ...
        '(%.1f to %.1f Angstrom)\n'],cleanup.number_pixels, ...
        numel(cleanup.affected_spectra),cleanup.minimum_wavelength, ...
        cleanup.maximum_wavelength)
else
    fprintf('No non-physical merged pixels were found\n')
end

end

function reject_file=find_merged_spectrum_rejections(fullwave,fullint)

reject_file=false(size(fullint,1),1);
keep_file=true(size(fullint,1),1);
maximum_rejection_iterations=50;

for rejection_iteration=1:maximum_rejection_iterations

    if sum(keep_file)<3
        break
    end

    medspec=median(fullint(keep_file,:),1,'omitnan');
    standdev=std(fullint(keep_file,:),0,1,'omitnan');

    bad_reference=~isfinite(medspec) | ~isfinite(standdev) | standdev==0;
    medspec(bad_reference)=1;
    standdev(bad_reference)=Inf;

    candidate_file=false(size(fullint,1),1);

    for file_index=1:size(fullint,1)

        if ~keep_file(file_index)
            continue
        end

        testspec=fullint(file_index,:);
        valid_test=isfinite(testspec);
        SS=valid_test & abs(testspec-medspec)>8*standdev;
        bad_fraction=sum(SS)/max(sum(valid_test),1);
        bad_wavelengths=fullwave(SS);
        if numel(bad_wavelengths)>1
            bad_span=max(bad_wavelengths)-min(bad_wavelengths);
        else
            bad_span=0;
        end
        valid_fraction=mean(valid_test);
        candidate_file(file_index)=valid_fraction<0.5 || ...
            (bad_fraction>0.1 && bad_span>1000);

    end

    new_reject=candidate_file & keep_file;

    if ~any(new_reject)
        fprintf('Merged-spectrum cleanup converged after %.0f iterations\n',rejection_iteration-1)
        break
    end

    reject_file=reject_file | new_reject;
    keep_file(new_reject)=false;

    fprintf('Merged-spectrum cleanup iteration %.0f flagged %.0f additional spectra\n',rejection_iteration,sum(new_reject))

    if rejection_iteration==maximum_rejection_iterations
        warning('Merged-spectrum cleanup reached the maximum number of rejection iterations before formal convergence.')
    end

end

end

function move_rejected_reduced_files(list,reject_indices,reduced_destination,processing_destination)

if ~isfolder(reduced_destination), mkdir(reduced_destination); end
if ~isfolder(processing_destination), mkdir(processing_destination); end

for nn=1:numel(reject_indices)

    list_index=reject_indices(nn);
    reduced_file_name=get_list_file_path(list(list_index));
    file=get_reduced_file_base_name(reduced_file_name);

    move_file_if_present(reduced_file_name,reduced_destination);
    prc_name=[file '_prc.mat'];
    move_file_if_present(fullfile(fileparts(processing_destination),prc_name),processing_destination);

end

end

function save_cleanup_report(report_file_name,list,cleanup_report,reject_file)

cleanup_manifest=struct();
cleanup_manifest.file_names={list.name};
cleanup_manifest.rejected=reject_file;
cleanup_manifest.rejected_file_names={list(reject_file).name};
cleanup_manifest.cleanup_report=cleanup_report;
cleanup_manifest.created_on=datestr(now);
save(report_file_name,'cleanup_manifest')

end

function delete_existing_chunk_files(chunk_directory)

existing_chunks=dir(fullfile(chunk_directory,'chunk_*.mat'));

for chunk_index=1:numel(existing_chunks)

    if ~existing_chunks(chunk_index).isdir
        delete(fullfile(existing_chunks(chunk_index).folder,existing_chunks(chunk_index).name))
    end

end

end

function chunk_size=get_effective_chunk_size(options,number_files)

if isfield(options,'chunked_post_reduction') && options.chunked_post_reduction==1
    chunk_size=options.post_reduction_chunk_size;
else
    chunk_size=number_files;
end

chunk_size=max(1,min(chunk_size,number_files));

end

function synthetic_vsini=get_synthetic_vsini(input_vsini)

%Choose vsini value for synthetic-spectrum and LSD line-mask generation
if isnumeric(input_vsini) && isscalar(input_vsini) && isfinite(input_vsini) && input_vsini>=0
    synthetic_vsini=input_vsini;
else
    synthetic_vsini=0;
    warning('Input vsini is missing or non-finite. Using vsini=0 km/s for synthetic-spectrum and LSD line-mask generation only.')
end

end

function primary_variable=load_primary_variable_from_mat_file(file_name)

%Load first variable from MAT file
loaded_file=load(file_name);
loaded_variable_names=fieldnames(loaded_file);
if isempty(loaded_variable_names), error(['MAT file contains no variables: ' file_name]); end
primary_variable=loaded_file.(loaded_variable_names{1});

end

function value=get_loaded_variable(loaded_file,preferred_variable_name)

%Get preferred or first variable
if isfield(loaded_file,preferred_variable_name)
    value=loaded_file.(preferred_variable_name);
else
    loaded_variable_names=fieldnames(loaded_file);
    if isempty(loaded_variable_names), error('Loaded MAT file contains no variables.'); end
    value=loaded_file.(loaded_variable_names{1});
end

end

function processing_frame=load_processing_frame(processing_directory,file)

%Load matching processing file
processing_file_name=fullfile(processing_directory,[file '_prc.mat']);
if ~isfile(processing_file_name), error(['Missing processing file: ' processing_file_name]); end
loaded_processing_file=load(processing_file_name);
expected_variable_name=[file '_prc'];
if isfield(loaded_processing_file,expected_variable_name)
    processing_frame=loaded_processing_file.(expected_variable_name);
else
    loaded_variable_names=fieldnames(loaded_processing_file);
    processing_frame=loaded_processing_file.(loaded_variable_names{1});
    warning(['Could not find expected variable ' expected_variable_name ' in ' processing_file_name '. Using first variable in file instead.']);
end

end

function [tmpwave,applied_bcorr,reduced_frame]=get_barycentric_corrected_wave(reduced_frame,objectname,file,ord)

%Apply or calculate barycentric correction
order_field_name=['order_' num2str(ord)];
if isfield(reduced_frame,'bcorr')
    stored_bcorr=reduced_frame.bcorr;
    reduced_frame.bcorr_was_updated=false;
else
    stored_bcorr=[];
end
if isempty(stored_bcorr)
    file_name=[file '.fit'];
    reduced_frame.bcorr_was_updated=true;
    if ~isfield(reduced_frame,'expdate'), error(['Missing expdate in ' file '.mat.']); end
    expdate=reduced_frame.expdate;
    if isfield(reduced_frame,'midtime_tt')
        midtime_tt=reduced_frame.midtime_tt;
    else
        if ~isfield(reduced_frame,'midtime_tt') && ~isfield(reduced_frame,'midtime_hrsp'), error(['Missing midtime_tt and midtime_hrsp in ' file '.mat.']); end
        midtime_utc=reduced_frame.midtime_hrsp;
        midtime_utc_dt=datetime([expdate(2:end-1) ' ' midtime_utc(2:end-1)],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        midtime_tt=midtime_utc_dt+seconds(69.184);%%%%%%%%%%%%%%%%update as needed
        midtime_tt=['"' datestr(midtime_tt,'HH:MM:SS.FFF') '"'];
        reduced_frame.midtime_utc=midtime_utc;
        reduced_frame.midtime_tt=midtime_tt;
    end
    stored_bcorr=str2double(get_barycentric_correction(objectname,expdate,midtime_tt,file_name));
    if isnan(stored_bcorr), stored_bcorr=0; end
    reduced_frame.bcorr=stored_bcorr;
end
applied_bcorr=stored_bcorr*-1;
if isempty(applied_bcorr)
    tmpwave=reduced_frame.(order_field_name)(:,1);
    applied_bcorr=0;
else
    tmpwave=shift_velocity(reduced_frame.(order_field_name)(:,1),applied_bcorr);
end

end

function save_reduced_frame_with_original_variable_name(file_name,variable_name,reduced_frame)

%Save updated reduced frame
saved_frame=struct();
saved_frame.(variable_name)=reduced_frame;
save(file_name,'-struct','saved_frame')

end

function move_file_if_present(source_file_name,destination_directory)

%Move file if it exists
if ~isfolder(destination_directory), mkdir(destination_directory); end
if isfile(source_file_name)
    [move_status,move_message]=movefile(source_file_name,destination_directory,'f');
    if ~move_status, warning(['Could not move ' source_file_name ': ' move_message]); end
end

end

function normalised_intensity=apply_continuum_vector_to_intensity_matrix(intensity_matrix,continuum_vector,intensity_field_name)

%Apply continuum vector to all spectra in an intensity matrix
continuum_vector=continuum_vector(:).';

if numel(continuum_vector)~=size(intensity_matrix,2)
    error(['Continuum fit length mismatch for ' intensity_field_name ': continuum has ' num2str(numel(continuum_vector)) ' points but intensity matrix has ' num2str(size(intensity_matrix,2)) ' wavelength pixels.'])
end

bad_continuum=~isfinite(continuum_vector) | continuum_vector==0;

if any(bad_continuum)
    warning(['Continuum fit for ' intensity_field_name ' contains non-finite or zero values. Replacing these values with 1.'])
    continuum_vector(bad_continuum)=1;
end

normalised_intensity=intensity_matrix./repmat(continuum_vector,size(intensity_matrix,1),1);

end

function regularised_lsd_options=get_regularised_lsd_options(options)

%Set regularised LSD options
regularised_lsd_options=struct();

if isfield(options,'extend_velocities') && options.extend_velocities==1
    regularised_lsd_options.velocity_limit=400;
else
    regularised_lsd_options.velocity_limit=200;
end

regularised_lsd_options.velocity_step=2.0;
regularised_lsd_options.regularisation_lambda=10;
regularised_lsd_options.regularisation_order=2;
regularised_lsd_options.normalisation='none';
regularised_lsd_options.profile_baseline='wing_median';

end

function error_values=make_error_vector(error_input,jd)

%Format errorbar values
if isscalar(error_input)
    error_values=error_input*ones(numel(jd),1);
else
    error_values=error_input;
end

end

function object_display_name=get_object_display_name(objectname)

%Format object label
if numel(objectname)>=4
    object_display_name=objectname(4:end);
else
    object_display_name=objectname;
end

end

function final_data_file_name=get_final_data_file_name(objectname,options)

%Choose final-data filename
if options.order_merge==1
    final_data_file_name=['final_data_' objectname '.mat'];
else
    final_data_file_name=['final_data_unmerged_' objectname '.mat'];
end

end

function reduced_files=list_reduced_spectrum_files(reduced_frames_directory)

%List reduced stellar files only
all_files=dir(fullfile(reduced_frames_directory,'J*.mat'));
keep_file=false(numel(all_files),1);
for file_index=1:numel(all_files)
    keep_file(file_index)=~isempty(regexp(all_files(file_index).name,'^J\d+\.mat$','once'));
end
reduced_files=all_files(keep_file);

end

function file_path=get_list_file_path(file_entry)

if isfield(file_entry,'folder') && ~isempty(file_entry.folder)
    file_path=fullfile(file_entry.folder,file_entry.name);
else
    file_path=file_entry.name;
end

end

function [finaldata,diagnostics]=separate_finaldata_diagnostics(finaldata)

diagnostics=struct('continuum_fitting',struct(),'rv_vsini',struct());
continuum_fields={'merged_suppnet','automatic_fit_method', ...
    'automatic_fit_summary','automatic_order_scales','template_fit_method', ...
    'template_fit_summary'};
rv_fields={'vbroad_interpretation','vsini_interpretation', ...
    'vsini_instrumental_fwhm', ...
    'vsini_instrumental_calibration_file','rv_vsini_moment_diagnostics'};

[finaldata,diagnostics.continuum_fitting]=move_structure_fields( ...
    finaldata,continuum_fields);
[finaldata,diagnostics.rv_vsini]=move_structure_fields(finaldata,rv_fields);

end

function [source,moved]=move_structure_fields(source,field_names)

moved=struct();
for field_index=1:numel(field_names)
    field_name=field_names{field_index};
    if isfield(source,field_name)
        moved.(field_name)=source.(field_name);
        source=rmfield(source,field_name);
    end
end

end

function diagnostics=merge_saved_diagnostics(diagnostics_file,new_diagnostics)

diagnostics=struct('continuum_fitting',struct(),'rv_vsini',struct());
if isfile(diagnostics_file)
    loaded=load(diagnostics_file,'diagnostics');
    if isfield(loaded,'diagnostics') && isstruct(loaded.diagnostics)
        diagnostics=loaded.diagnostics;
    end
end

sections={'continuum_fitting','rv_vsini'};
for section_index=1:numel(sections)
    section=sections{section_index};
    if ~isfield(diagnostics,section), diagnostics.(section)=struct(); end
    if ~isfield(new_diagnostics,section), continue; end
    fields=fieldnames(new_diagnostics.(section));
    for field_index=1:numel(fields)
        field_name=fields{field_index};
        diagnostics.(section).(field_name)= ...
            new_diagnostics.(section).(field_name);
    end
end

end

function sorted_files=sort_reduced_spectrum_files(reduced_files)

observation_numbers=nan(numel(reduced_files),1);

for file_index=1:numel(reduced_files)

    file_base=get_reduced_file_base_name(reduced_files(file_index).name);
    observation_numbers(file_index)=str2double(file_base(2:end));

end

[~,sort_index]=sort(observation_numbers);
sorted_files=reduced_files(sort_index);

end

function file_base=get_reduced_file_base_name(file_name)

if isstring(file_name)
    file_name=char(file_name);
end

[~,file_name,file_extension]=fileparts(file_name);
file_name=[file_name file_extension];

file_match=regexp(file_name,'^(J\d+)','tokens','once');

if isempty(file_match)
    error(['Could not extract J-number from file name: ' file_name])
end

file_base=file_match{1};

end

function output=rounding(input,tonearest)

%Round input to nearest requested value
output=round(input/tonearest)*tonearest;

end
