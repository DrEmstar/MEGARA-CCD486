function [summary,line_measurements,order_measurements]=calculate_resolving_power( ...
    thar_images,flat_info_file,blue_data_chop_value,thar_info_file,options)
%Calculate resolving power from isolated ThAr lines.
%
%This is an optional diagnostic and is not called during a normal run.
%
%Example:
% [summary,lines,orders]=calculate_resolving_power('J9236023.fit');

%Put the MEGARA 2.1 reduction functions first on the path.
this_file=mfilename('fullpath');
version_root=fileparts(fileparts(this_file));
addpath(genpath(fullfile(version_root,'Reduction_code')),'-begin')

%Check the input files.
if nargin<1 || isempty(thar_images)
    error('Supply at least one CCD486 ThAr FITS file.')
end
if ischar(thar_images) || (isstring(thar_images) && isscalar(thar_images))
    thar_images={char(thar_images)};
elseif isstring(thar_images)
    thar_images=cellstr(thar_images(:));
elseif iscell(thar_images)
    thar_images=cellfun(@char,thar_images(:),'UniformOutput',false);
else
    error('thar_images must contain one or more filenames.')
end
for file_index=1:numel(thar_images)
    if ~isfile(thar_images{file_index})
        error(['ThAr FITS file not found: ' thar_images{file_index}])
    end
end

%Use the reduction files beside the first ThAr unless otherwise supplied.
first_thar_file=get_full_filename(thar_images{1});
first_thar_folder=fileparts(first_thar_file);
if nargin<3 || isempty(blue_data_chop_value)
    blue_data_chop_value=300;
end
if nargin<2 || isempty(flat_info_file)
    flat_info_file=fullfile(first_thar_folder, ...
        ['flat_info_blue_' num2str(blue_data_chop_value) '.mat']);
end
if nargin<4 || isempty(thar_info_file)
    thar_info_file=fullfile(first_thar_folder, ...
        ['ThAr_info_blue_' num2str(blue_data_chop_value) '.mat']);
end
if nargin<5 || isempty(options)
    options=struct();
end
options=set_default_options(options,first_thar_folder);

if ~isfile(flat_info_file)
    error(['Flat information file not found: ' flat_info_file])
end
if ~isfile(thar_info_file)
    error(['ThAr information file not found: ' thar_info_file])
end

%Load the order traces, extraction profiles and wavelength solutions.
flat_information=load(flat_info_file);
if ~isfield(flat_information,'allorders') || ~isfield(flat_information,'flatdata')
    error([flat_info_file ' must contain allorders and flatdata.'])
end
allorders=flat_information.allorders;
flatdata=flat_information.flatdata;

thar_information=load(thar_info_file);
if ~isfield(thar_information,'jd_ths') || ~isfield(thar_information,'fwave')
    error([thar_info_file ' must contain jd_ths and fwave.'])
end

%Measure each requested exposure.
line_tables=cell(numel(thar_images),1);
frame_information=cell(numel(thar_images),1);
for file_index=1:numel(thar_images)
    [line_tables{file_index},frame_information{file_index}]=measure_one_frame( ...
        thar_images{file_index},blue_data_chop_value,allorders,flatdata, ...
        thar_information,options);
end
line_measurements=vertcat(line_tables{:});

%Remove outliers independently in each order.
line_measurements.AcceptedBeforeRobustClip=line_measurements.Accepted;
absolute_orders=unique(line_measurements.Order(isfinite(line_measurements.Order)));
for order_index=1:numel(absolute_orders)
    order_number=absolute_orders(order_index);
    indices=find(line_measurements.Order==order_number & line_measurements.Accepted);
    if numel(indices)<options.minimum_lines_per_order
        continue
    end
    values=line_measurements.ResolvingPower(indices);
    centre=median(values,'omitnan');
    spread=median(abs(values-centre),'omitnan');
    if isfinite(spread) && spread>0
        bad=abs(values-centre)>options.robust_clip_mad*1.4826*spread;
        bad_indices=indices(bad);
        line_measurements.Accepted(bad_indices)=false;
        line_measurements.RejectionReason(bad_indices)="robust order outlier";
    end
end

%Make the order and overall summaries.
order_measurements=make_order_summary(line_measurements, ...
    options.minimum_lines_per_order);
accepted_values=line_measurements.ResolvingPower(line_measurements.Accepted);
if isempty(accepted_values)
    error('No ThAr lines passed the resolving-power checks.')
end

summary=struct();
summary.version='MEGARA_2.0_RESOLVING_POWER_1';
summary.created_utc=char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
summary.thar_images=thar_images;
summary.flat_info_file=flat_info_file;
summary.thar_info_file=thar_info_file;
summary.blue_data_chop_value=blue_data_chop_value;
summary.options=options;
summary.frame_information=[frame_information{:}];
summary.reference_files=unique(string({summary.frame_information.reference_file}));
summary.number_fitted_lines=height(line_measurements);
summary.number_accepted_lines=nnz(line_measurements.Accepted);
summary.number_orders=nnz(order_measurements.NumberAcceptedLines>= ...
    options.minimum_lines_per_order);
summary.median_resolving_power=median(accepted_values,'omitnan');
summary.mean_resolving_power=mean(accepted_values,'omitnan');
summary.resolving_power_16_84=prctile(accepted_values,[16 84]);
summary.median_velocity_fwhm_kms=299792.458/summary.median_resolving_power;

fprintf('\nMEGARA 2.1 CCD486 resolving-power diagnostic\n')
fprintf('Accepted %d/%d fitted lines in %d orders.\n', ...
    summary.number_accepted_lines,summary.number_fitted_lines,summary.number_orders)
fprintf('Median R = %.0f (16th-84th percentile: %.0f-%.0f).\n', ...
    summary.median_resolving_power,summary.resolving_power_16_84(1), ...
    summary.resolving_power_16_84(2))
fprintf('Equivalent median FWHM = %.2f km/s.\n', ...
    summary.median_velocity_fwhm_kms)

if options.make_plots
    make_diagnostic_plots(line_measurements,order_measurements,summary)
end
if options.save_results
    save(options.output_file,'summary','line_measurements','order_measurements')
    fprintf('Saved resolving-power results: %s\n',options.output_file)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [measurements,frame_information]=measure_one_frame(filename,blue_chop, ...
    allorders,flatdata,thar_information,options)

filename=get_full_filename(filename);

%Read the CCD486 frame in its established orientation.
data=hercules_fitsread(filename);
data=blue_data_chop(data,blue_chop);
extracted=extract_all_orders_no_background(allorders,data);
extracted=apply_flat_profile_weights(extracted,flatdata);
[thardata,~]=redo_order_numbering_simple(extracted,allorders);

%Get the observation time and historical line catalogue.
header=get_header(filename)';
header=get_required_headers_from_header(header);
header=check_hercules_headers(header);
if ~isfield(header,'MJD_OBS') || ~isfinite(header.MJD_OBS)
    error(['No usable MJD_OBS/JD was recovered from ' filename])
end
[reference,reference_file]=load_reference_file(filename,header.MJD_OBS);

%Use the wavelength solution nearest this exposure.
[~,wavelength_solution_index]=min(abs(thar_information.jd_ths-header.MJD_OBS));
wavelength_field=['th' num2str(wavelength_solution_index)];
if ~isfield(thar_information.fwave,wavelength_field)
    error(['Saved wavelength solution ' wavelength_field ' is missing.'])
end
wave=thar_information.fwave.(wavelength_field);

%Fit the known ThAr lines using the normal wavelength-calibration fitter.
[thfitinfo,~,move_file_flag]=find_th_lines(thardata,reference.order, ...
    reference.ul,reference.uc,reference.ur,0);
if move_file_flag || isempty(thfitinfo)
    error(['Too few usable ThAr lines were found in ' filename])
end

number_lines=size(thfitinfo,1);
file_column=repmat(string(filename),number_lines,1);
cache_index_column=repmat(wavelength_solution_index,number_lines,1);
reference_index=round(thfitinfo(:,1));
order_column=double(reference.order(reference_index));
catalogue_wavelength=double(reference.air(reference_index));
centroid_pixel=thfitinfo(:,4);
fwhm_pixel=thfitinfo(:,5);
fit_resnorm=thfitinfo(:,2);
fitted_height=thfitinfo(:,3);
measured_wavelength=nan(number_lines,1);
dispersion=nan(number_lines,1);
peak_counts=nan(number_lines,1);
nearest_separation=nan(number_lines,1);

%Convert each fitted pixel width to wavelength.
for line_index=1:number_lines
    order_number=order_column(line_index);
    order_field=['order_' num2str(order_number)];
    if ~isfield(thardata,order_field) || ~isfield(wave,order_field)
        continue
    end
    x=double(thardata.(order_field).xax(:));
    flux=double(thardata.(order_field).summed_data(:));
    wavelength=double(wave.(order_field)(:));
    if numel(x)~=numel(wavelength) || numel(x)~=numel(flux)
        continue
    end
    measured_wavelength(line_index)=interp1(x,wavelength, ...
        centroid_pixel(line_index),'linear',NaN);
    local_dispersion=gradient(wavelength)./gradient(x);
    dispersion(line_index)=interp1(x,local_dispersion, ...
        centroid_pixel(line_index),'linear',NaN);
    local_pixels=abs(x-centroid_pixel(line_index))<=1;
    if any(local_pixels)
        peak_counts(line_index)=max(flux(local_pixels),[],'omitnan');
    end
    same_order=find(reference.order==order_number);
    other_centres=abs(reference.uc(same_order)- ...
        reference.uc(reference_index(line_index)));
    other_centres=other_centres(other_centres>0);
    if ~isempty(other_centres)
        nearest_separation(line_index)=min(other_centres);
    end
end

fwhm_wavelength=abs(fwhm_pixel.*dispersion);
resolving_power=measured_wavelength./fwhm_wavelength;
velocity_fwhm=299792.458./resolving_power;
accepted=true(number_lines,1);
reason=strings(number_lines,1);

invalid=~isfinite(resolving_power) | ~isfinite(measured_wavelength) | dispersion==0;
[accepted,reason]=reject_lines(accepted,reason,invalid,"invalid wavelength conversion");
[accepted,reason]=reject_lines(accepted,reason, ...
    fit_resnorm>options.maximum_fit_resnorm,"poor Gaussian fit");
[accepted,reason]=reject_lines(accepted,reason, ...
    fwhm_pixel<options.minimum_fwhm_pixels,"FWHM below limit");
[accepted,reason]=reject_lines(accepted,reason, ...
    fwhm_pixel>options.maximum_fwhm_pixels,"FWHM above limit");
[accepted,reason]=reject_lines(accepted,reason, ...
    resolving_power<options.minimum_resolving_power,"R below limit");
[accepted,reason]=reject_lines(accepted,reason, ...
    resolving_power>options.maximum_resolving_power,"R above limit");
blend=nearest_separation<options.minimum_isolation_fwhm.*fwhm_pixel;
[accepted,reason]=reject_lines(accepted,reason,blend,"nearby catalogue line");

measurements=table(file_column,cache_index_column,reference_index,order_column, ...
    catalogue_wavelength,measured_wavelength,centroid_pixel,fwhm_pixel, ...
    dispersion,fwhm_wavelength,resolving_power,velocity_fwhm,peak_counts, ...
    fitted_height,fit_resnorm,nearest_separation,accepted,reason, ...
    'VariableNames',{'File','WavelengthCacheIndex','ReferenceIndex','Order', ...
    'CatalogueWavelength_A','MeasuredWavelength_A','CentroidPixel','FWHM_Pixel', ...
    'Dispersion_A_PerPixel','FWHM_A','ResolvingPower','VelocityFWHM_kms', ...
    'PeakExtractedCounts','FittedNormalisedHeight','FitResidualNorm', ...
    'NearestReferenceSeparation_Pixel','Accepted','RejectionReason'});

frame_information=struct();
frame_information.file=filename;
frame_information.header_jd=header.MJD_OBS;
frame_information.wavelength_cache_index=wavelength_solution_index;
frame_information.wavelength_cache_jd= ...
    thar_information.jd_ths(wavelength_solution_index);
frame_information.reference_file=reference_file;
frame_information.number_fitted_lines=number_lines;
frame_information.number_initially_accepted_lines=nnz(accepted);

end

function extracted=apply_flat_profile_weights(extracted,flatdata)

%Use the same flat-profile weighting as the normal ThAr reduction.
for order_index=1:flatdata.numords
    field=['order_' num2str(order_index)];
    if ~isfield(extracted,field) || ~isfield(flatdata,field)
        continue
    end
    indices=flatdata.(field).extraction_inds;
    if isfield(flatdata.(field),'profile_weights')
        weights=flatdata.(field).profile_weights;
    else
        weights=zeros(size(indices));
        weights(indices)=1;
        weights=weights/sum(weights);
    end
    extracted.(field).summed_data=sum( ...
        extracted.(field).data(:,indices).*weights(indices),2);
end

end

function [reference,reference_file]=load_reference_file(filename,jd)

%Select the historical catalogue used by the normal MEGARA 2.1 reduction.
[~,short_name]=fileparts(filename);
if numel(short_name)<5
    error(['Could not determine the date code from ' filename])
end
thardate=str2double(short_name(2:5));
if ~isfinite(thardate)
    error(['Could not determine the date code from ' filename])
end
reference_name=select_reference_name(thardate,jd);
reference_file=which(reference_name);
if isempty(reference_file)
    error(['Could not find the historical ThAr file ' reference_name])
end
reference=load(reference_file);
required_fields={'order','air','ul','uc','ur'};
for field_index=1:numel(required_fields)
    if ~isfield(reference,required_fields{field_index})
        error([reference_file ' is missing ' required_fields{field_index} '.'])
    end
end

end

function reference_name=select_reference_name(thardate,jd)

%Keep the established MEGARA 2.1 date ranges.
if jd<60000
    if thardate<=4463
        reference_name='thardat_g.mat';
    elseif thardate<=4781
        reference_name='thardat_f.mat';
    elseif thardate<=5454
        reference_name='thardat_d.mat';
    elseif thardate<=5702
        reference_name='thardat_e.mat';
    elseif thardate<=6337
        reference_name='thardat_c.mat';
    elseif thardate<6791
        reference_name='thardat_b.mat';
    elseif thardate<7083
        reference_name='thardat_a.mat';
    elseif thardate<7178
        reference_name='thardat_i.mat';
    elseif thardate<7821
        reference_name='thardat_a.mat';
    elseif thardate<7968
        reference_name='thardat_h.mat';
    elseif thardate<8120
        reference_name='thardat_a.mat';
    elseif thardate<8144
        reference_name='thardat_j.mat';
    elseif thardate<8146
        reference_name='thardat_k.mat';
    elseif thardate<8303
        reference_name='thardat_j.mat';
    elseif thardate<9461
        reference_name='thardat_l.mat';
    elseif thardate<9650
        reference_name='thardat_mm.mat';
    else
        reference_name='thardat_n.mat';
    end
elseif thardate<889
    reference_name='thardat_o.mat';
else
    reference_name='thardat_p.mat';
end

end

function [accepted,reason]=reject_lines(accepted,reason,reject_mask,message)

new_rejections=accepted & reject_mask;
accepted(new_rejections)=false;
reason(new_rejections)=message;

end

function order_table=make_order_summary(lines,minimum_lines)

orders=unique(lines.Order(isfinite(lines.Order)));
number_orders=numel(orders);
number_fitted=zeros(number_orders,1);
number_accepted=zeros(number_orders,1);
median_wavelength=nan(number_orders,1);
median_fwhm_pixels=nan(number_orders,1);
median_r=nan(number_orders,1);
r16=nan(number_orders,1);
r84=nan(number_orders,1);
for index=1:number_orders
    in_order=lines.Order==orders(index);
    accepted=in_order & lines.Accepted;
    number_fitted(index)=nnz(in_order);
    number_accepted(index)=nnz(accepted);
    if number_accepted(index)>=minimum_lines
        median_wavelength(index)=median(lines.MeasuredWavelength_A(accepted),'omitnan');
        median_fwhm_pixels(index)=median(lines.FWHM_Pixel(accepted),'omitnan');
        median_r(index)=median(lines.ResolvingPower(accepted),'omitnan');
        percentiles=prctile(lines.ResolvingPower(accepted),[16 84]);
        r16(index)=percentiles(1);
        r84(index)=percentiles(2);
    end
end
order_table=table(orders,number_fitted,number_accepted,median_wavelength, ...
    median_fwhm_pixels,median_r,r16,r84,'VariableNames', ...
    {'Order','NumberFittedLines','NumberAcceptedLines','MedianWavelength_A', ...
    'MedianFWHM_Pixel','MedianResolvingPower','ResolvingPower16', ...
    'ResolvingPower84'});

end

function make_diagnostic_plots(lines,orders,summary)

accepted=lines.Accepted;
figure('Color','w','Name','MEGARA 2.1 resolving power', ...
    'Position',[100 100 1400 850])
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

nexttile
scatter(lines.MeasuredWavelength_A(accepted),lines.ResolvingPower(accepted), ...
    12,lines.Order(accepted),'filled','MarkerFaceAlpha',0.45)
yline(summary.median_resolving_power,'k--','Median')
xlabel('Wavelength (Angstrom)')
ylabel('Resolving power, R')
title('Accepted isolated ThAr lines')
colour_bar=colorbar;
colour_bar.Label.String='Absolute order';
grid on

nexttile
valid=orders.NumberAcceptedLines>=summary.options.minimum_lines_per_order;
errorbar(orders.Order(valid),orders.MedianResolvingPower(valid), ...
    orders.MedianResolvingPower(valid)-orders.ResolvingPower16(valid), ...
    orders.ResolvingPower84(valid)-orders.MedianResolvingPower(valid), ...
    'o','MarkerFaceColor',[0.2 0.45 0.8])
yline(summary.median_resolving_power,'k--')
xlabel('Absolute order')
ylabel('Median resolving power, R')
title('Order-by-order median')
grid on

nexttile
histogram(lines.ResolvingPower(accepted),40)
xline(summary.median_resolving_power,'k--','Median')
xlabel('Resolving power, R')
ylabel('Number of lines')
title('Accepted-line distribution')
grid on

nexttile
scatter(lines.MeasuredWavelength_A(accepted),lines.FWHM_Pixel(accepted), ...
    12,lines.Order(accepted),'filled','MarkerFaceAlpha',0.45)
xlabel('Wavelength (Angstrom)')
ylabel('Fitted FWHM (pixels)')
title('Detector sampling')
grid on

sgtitle(sprintf('MEGARA 2.1 CCD486: median R = %.0f from %d lines', ...
    summary.median_resolving_power,summary.number_accepted_lines), ...
    'FontWeight','bold')

end

function options=set_default_options(options,output_folder)

if ~isstruct(options)
    error('options must be a structure.')
end
defaults.make_plots=true;
defaults.save_results=false;
output_time=char(datetime('now','Format','yyyyMMdd_HHmmss'));
defaults.output_file=fullfile(output_folder, ...
    ['resolving_power_' output_time '.mat']);
defaults.maximum_fit_resnorm=5;
defaults.minimum_fwhm_pixels=1.6;
defaults.maximum_fwhm_pixels=8.5;
defaults.minimum_resolving_power=20000;
defaults.maximum_resolving_power=200000;
defaults.minimum_isolation_fwhm=2.5;
defaults.robust_clip_mad=4;
defaults.minimum_lines_per_order=3;
field_names=fieldnames(defaults);
for index=1:numel(field_names)
    field=field_names{index};
    if ~isfield(options,field) || isempty(options.(field))
        options.(field)=defaults.(field);
    end
end

end

function filename=get_full_filename(filename)

if isfile(filename)
    file_information=dir(filename);
    filename=fullfile(file_information(1).folder,file_information(1).name);
    return
end
path_match=which(filename);
if ~isempty(path_match)
    filename=path_match;
end

end
