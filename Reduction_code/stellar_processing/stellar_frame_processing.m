function stellar_frame_processing(stellarimgname,info,allorders,flatdata,fwave1,fwave2,jd1,jd2,blue_data_chop_value,flatdata_medfilt,backdata_medfilt,cosmic_plotting,jdexp)
%STELLAR image processing

%Revised Apr2026- median taken only from the trimmed, well-illuminated part of the order.

%Initialise outputs
reduced_spectrum_data=struct();
processing_data=struct();

%Get file root
sname=stellarimgname(1:end-4);

%Read Julian date
if isempty(jdexp)
    jdexp=get_julian_date_from_info(info,sname);
end
if isempty(jdexp) || ~isfinite(jdexp)
    disp('couldn''t get Julian date, check filename and')
    disp('info file and try again')
    return
end

%Read and prepare stellar image
data=raw_hercules_fitsread(stellarimgname);
data=blue_data_chop(data,blue_data_chop_value);

%Fit stellar background
backfit=fit_background(data,allorders);

%Register science extraction to this exposure's monthly flat profile.
[stardata,backdata,trace_registration]=extract_registered_stellar_orders( ...
    allorders,data,backfit,flatdata,cosmic_plotting);
processing_data.trace_registration=trace_registration;
reduced_spectrum_data.trace_registration=trace_registration;

%Apply flat-profile weighted extraction
[stardata,backdata]=apply_flat_profile_weighted_extraction(stardata,backdata,flatdata,backdata_medfilt);

%apply simple program to change to absolute order numbers
[stardata,m]=redo_order_numbering_simple(stardata,allorders);

%make wavelength axes for stellar obs by linear interpolation of ThAr axes
wavelength=linear_interpolate_wavelengths1D(fwave1,fwave2,jd1,jd2,jdexp,m);

%Set trimming and normalisation thresholds
trim_threshold=0.1;
norm_threshold=0.5;
flat_quality_floor=0.2;
flat_quality_full=0.5;

%Flat-field and normalise each order
for order_index=1:numel(m)
    original_order_number=order_index;
    absolute_order_number=m(order_index);
    flat_order_field_name=['order_' num2str(original_order_number)];
    stellar_order_field_name=['order_' num2str(absolute_order_number)];
    if ~isfield(flatdata,flat_order_field_name), warning(['Skipping order ' num2str(absolute_order_number) ' because flatdata.' flat_order_field_name ' is missing.']); continue; end
    if ~isfield(stardata,stellar_order_field_name), warning(['Skipping order ' num2str(absolute_order_number) ' because stardata.' stellar_order_field_name ' is missing.']); continue; end
    if ~isfield(wavelength,stellar_order_field_name), warning(['Skipping order ' num2str(absolute_order_number) ' because wavelength.' stellar_order_field_name ' is missing.']); continue; end
    fl=flatdata.(flat_order_field_name).sumdata;
    flb=fl;
    fl=medfilt1(fl,flatdata_medfilt); %median filter the flat-field data
    fl=fl(:);
    flat_normalisation=median(fl(isfinite(fl) & fl>0));
    if ~isfinite(flat_normalisation) || flat_normalisation==0, flat_normalisation=1; end
    fl=fl/flat_normalisation;
    bad_flat_pixels=~isfinite(fl) | fl<=0;
    fl(bad_flat_pixels)=NaN;
    flat_quality=zeros(size(fl));
    valid_flat=~bad_flat_pixels;
    flat_quality(valid_flat)=(fl(valid_flat)-flat_quality_floor)./ ...
        (flat_quality_full-flat_quality_floor);
    flat_quality=min(max(flat_quality,0),1);
    xxx=stardata.(stellar_order_field_name).summed_data_b(:);
    xxx(~isfinite(xxx))=NaN;
    source_counts=stardata.(stellar_order_field_name).summed_data(:);
    background_counts=backdata.(flat_order_field_name).summed_smoothed_data(:);

    % trim based on actual threshold crossings in flat profile
    number_flat_pixels=numel(fl);
    midpoint=floor(number_flat_pixels/2);
    left_inds=find(fl(1:midpoint)>trim_threshold);
    right_inds=find(fl(midpoint+1:end)>trim_threshold);
    if isempty(left_inds)
        ind1=1;
    else
        ind1=left_inds(1);
    end
    if isempty(right_inds)
        ind2=number_flat_pixels;
    else
        ind2=midpoint+right_inds(end);
    end
    if ind2<=ind1
        ind1=1;
        ind2=number_flat_pixels;
    end
    trim_inds=ind1:ind2;

    % robust per-order normalisation using only well-illuminated part
    norm_mask=false(size(xxx));
    norm_mask(trim_inds)=true;
    norm_mask=norm_mask & (fl>norm_threshold) & isfinite(xxx) & isfinite(fl);
    if nnz(norm_mask)>20
        xmed=median(xxx(norm_mask));
    else
        xmed=median(xxx(trim_inds),'omitnan');
    end
    if ~isfinite(xmed) || xmed==0
        xmed=1;
    end
    xxx=xxx/xmed;
    xxxf=xxx(:)./fl(:);
    xxxf(flat_quality==0)=NaN;
    normalised_uncertainty=sqrt(max(source_counts,0)+ ...
        abs(background_counts))./(abs(xmed).*fl);
    extraction_quality=calculate_extraction_quality( ...
        normalised_uncertainty);
    xxxf(extraction_quality==0)=NaN;

    %Store reduced order
    reduced_order=cat(2,wavelength.(stellar_order_field_name)(:),xxxf(:));
    reduced_spectrum_data.(stellar_order_field_name)=reduced_order(ind1:ind2,:);

    %Store processing products
    processing_data.order_cutpixels.(stellar_order_field_name).ind1=ind1;
    processing_data.order_cutpixels.(stellar_order_field_name).ind2=ind2;
    processing_data.(['order_bkgd_' num2str(absolute_order_number)])=backdata.(flat_order_field_name).summed_data(:);
    processing_data.(['order_bkgd_smoothed_' num2str(absolute_order_number)])=backdata.(flat_order_field_name).summed_smoothed_data(:);
    processing_data.(['order_flat_' num2str(absolute_order_number)])=flb(:);
    processing_data.(['order_flat_smoothed_' num2str(absolute_order_number)])=fl(:);
    processing_data.(['order_flat_quality_' num2str(absolute_order_number)])=flat_quality(:);
    processing_data.(['order_normalised_uncertainty_' num2str(absolute_order_number)])=normalised_uncertainty(:);
    processing_data.(['order_extraction_quality_' num2str(absolute_order_number)])=extraction_quality(:);
    processing_data.(['order_nff_' num2str(absolute_order_number)])=stardata.(stellar_order_field_name).summed_data(:);
    processing_data.notes=['Reduced spectrum is: ($file_prc.order_nff_ - $file_prc.order_bkgd_smoothed_) / ($file_prc.order_flat_smoothed_). Invalid flat pixels are NaN; order_flat_quality_ tapers from 0 at flat response ' num2str(flat_quality_floor) ' to 1 at ' num2str(flat_quality_full) '. background smoothing is at medfilt =' num2str(backdata_medfilt) ', flat smoothing is at ' num2str(flatdata_medfilt)];
end

%Store Julian date
reduced_spectrum_data.jd=jdexp;

%Read stellar headers
header=get_header([sname '.fit']);
headers=get_required_headers_from_header(header);
starname=get_optional_header_value(headers,'OBJECT','');
expdate=get_exposure_date(headers);
starttime=get_optional_header_value(headers,'REC_STRT','');

%calculating signal-to-noise
exp_time=get_optional_header_value(headers,'EXPTIME',NaN);
signal_to_noise=calculate_stellar_signal_to_noise(stardata,exp_time);
reduced_spectrum_data.signal_to_noise=signal_to_noise;

%Adding the mid exposure time
midexp=get_optional_header_value(headers,'HERCFWMT',NaN);
if isempty(midexp) || isnan(midexp)
    warning('Flux-weighted mid-time unavailable for %s, using exposure mid-time instead',sname)
    midexp=exp_time/2;
end

%Calculate mid-times
[midtime_utc,midtime_tt]=calculate_mid_exposure_times(expdate,starttime,midexp,sname);

%saves some headers for adding to the '.mats'
fibre=get_optional_header_value(headers,'HERCFIB',NaN);     %Hercules Fibre 1/2/3
counts=get_optional_header_value(headers,'HERCFTC',NaN);     %Hercules flux meter total counts
RA_j2000=get_preferred_coordinate(headers,'RA','RA_2000','RA information missing'); %RA in j2000
DEC_j2000=get_preferred_coordinate(headers,'DEC','DEC_2000','DEC information missing'); %DEC in j2000

%Store metadata
reduced_spectrum_data.starname=starname;
reduced_spectrum_data.expdate=expdate;
reduced_spectrum_data.midtime_utc=midtime_utc;
reduced_spectrum_data.midtime_tt=midtime_tt;
reduced_spectrum_data.sname=sname;
reduced_spectrum_data.blue_data_chop=blue_data_chop_value;
reduced_spectrum_data.counts=counts;
reduced_spectrum_data.fibre=fibre;
reduced_spectrum_data.exptime=exp_time;
reduced_spectrum_data.RA_j2000=RA_j2000;
reduced_spectrum_data.DEC_j2000=DEC_j2000;

%put barycentric correction in .fits file
bcorr=str2double(get_barycentric_correction(starname,expdate,midtime_tt,stellarimgname));
reduced_spectrum_data.bcorr=bcorr;

%Save reduced and processing files
fprintf('Reduced %s.mat and %s.mat\n\n',sname,[sname '_prc'])
save_legacy_named_variable(sname,sname,reduced_spectrum_data)
save_legacy_named_variable([sname '_prc'],[sname '_prc'],processing_data)

% Get list of variable names in current function workspace
% Clear them all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extraction_quality=calculate_extraction_quality(normalised_uncertainty)

%Taper uncertain background-subtracted pixels to zero merge weight
full_quality_limit=0.10;
zero_quality_limit=0.20;
extraction_quality=(zero_quality_limit-normalised_uncertainty)./ ...
    (zero_quality_limit-full_quality_limit);
extraction_quality=min(max(extraction_quality,0),1);
extraction_quality(~isfinite(normalised_uncertainty))=0;

end

function jdexp=get_julian_date_from_info(info,sname)

%Read Julian date from info
if isfield(info,sname) && isfield(info.(sname),'MJD_OBS')
    jdexp=info.(sname).MJD_OBS;
else
    jdexp=[];
end

end

function [stardata,backdata]=apply_flat_profile_weighted_extraction(stardata,backdata,flatdata,backdata_medfilt)

%Apply weighted extraction order by order
if ~isfield(flatdata,'numords'), error('flatdata does not contain numords.'); end
for order_number=1:flatdata.numords
    order_field_name=['order_' num2str(order_number)];
    if ~isfield(flatdata,order_field_name), warning(['Skipping order ' num2str(order_number) ' because flatdata.' order_field_name ' is missing.']); continue; end
    if ~isfield(stardata,order_field_name), warning(['Skipping order ' num2str(order_number) ' because stardata.' order_field_name ' is missing.']); continue; end
    if ~isfield(backdata,order_field_name), warning(['Skipping order ' num2str(order_number) ' because backdata.' order_field_name ' is missing.']); continue; end
    ind=flatdata.(order_field_name).extraction_inds;
    profile=get_profile_weights(flatdata.(order_field_name),ind);
    bdata=backdata.(order_field_name).data;
    odata=stardata.(order_field_name).data;
    if numel(ind)~=size(odata,2) || numel(ind)~=size(bdata,2), error(['Extraction index length does not match extracted data width for ' order_field_name '.']); end
    backdata.(order_field_name).summed_data=sum(bdata(:,ind).*profile(ind),2);
    backdata.(order_field_name).summed_smoothed_data=mean_smoothing(backdata.(order_field_name).summed_data,backdata_medfilt);
    stardata.(order_field_name).summed_data=sum(odata(:,ind).*profile(ind),2);
    stardata.(order_field_name).summed_data_b=stardata.(order_field_name).summed_data-backdata.(order_field_name).summed_smoothed_data;
end

end

function profile=get_profile_weights(flat_order_data,ind)

%Read or build profile weights
if isfield(flat_order_data,'profile_weights')
    profile=flat_order_data.profile_weights;
else
    profile=zeros(size(ind));
    profile(ind)=1;
    profile_sum=sum(profile);
    if profile_sum==0 || ~isfinite(profile_sum), profile=ones(size(ind)); profile_sum=sum(profile); end
    profile=profile./profile_sum;
end
profile=profile(:)';
if numel(profile)~=numel(ind), error('profile_weights and extraction_inds must have the same length.'); end

end

function value=get_optional_header_value(headers,field_name,default_value)

%Read optional header value
if isfield(headers,field_name) && ~isempty(headers.(field_name))
    value=headers.(field_name);
else
    value=default_value;
end

end

function expdate=get_exposure_date(headers)

%Read exposure date
expdate=get_optional_header_value(headers,'DATE_OBS','');  % OLD CCD SOFTWARE
if numel(expdate)>=12 && expdate(12)=='T'
    expdate=get_optional_header_value(headers,'DATE',expdate);% NEW CCD SOFTWARE
end
if isempty(expdate)
    expdate=get_optional_header_value(headers,'DATE','');
end
if isempty(expdate), error('No DATE_OBS or DATE header found.'); end

end

function signal_to_noise=calculate_stellar_signal_to_noise(stardata,exp_time)

%Calculate signal-to-noise from preferred order
signal_to_noise=NaN;
if isfield(stardata,'order_100')
    stellar_order=stardata.order_100;
elseif isfield(stardata,'order_101')
    stellar_order=stardata.order_101;
else
    order_fields=fieldnames(stardata);
    order_fields=order_fields(startsWith(order_fields,'order_'));
    if isempty(order_fields), warning('No stellar order available for signal-to-noise calculation.'); return; end
    stellar_order=stardata.(order_fields{1});
end
if ~isfinite(exp_time), warning('Exposure time unavailable for signal-to-noise calculation.'); return; end
signal_to_noise=compute_signal_to_noise_spectrum(stellar_order,exp_time);

end

function [midtime_utc,midtime_tt]=calculate_mid_exposure_times(expdate,starttime,midexp,sname)

%Clean date and time strings
clean_expdate=clean_header_string(expdate);
clean_starttime=clean_header_string(starttime);
if isempty(clean_expdate) || isempty(clean_starttime), warning('Could not calculate mid-exposure time for %s.',sname); midtime_utc='""'; midtime_tt='""'; return; end

%Parse date and time
if contains(clean_expdate,'T')
    date_parts=split(clean_expdate,'T');
    clean_expdate=strtrim(date_parts{1});
end
try
    exposure_date=datetime(clean_expdate,'InputFormat','yyyy-MM-dd');
catch
    try
        exposure_date=datetime(clean_expdate,'InputFormat','dd/MM/yyyy');
    catch
        try
            exposure_date=datetime(clean_expdate(1:10),'InputFormat','yyyy-MM-dd');
        catch
            error('Could not interpret exposure date for %s: %s',sname,clean_expdate)
        end
    end
end
start_time_duration=normalise_header_time_duration(clean_starttime);
midtime=exposure_date+start_time_duration+seconds(midexp);

%putting in Terrestial Time instead of UTC for hrsp_barycorr
midtime_tt_datetime=midtime+seconds(69.184);
midtime_utc=['"' datestr(midtime,'HH:MM:SS.FFF') '"'];
midtime_tt=['"' datestr(midtime_tt_datetime,'HH:MM:SS.FFF') '"'];

end

function normalised_time_duration=normalise_header_time_duration(input_time)

input_time=clean_header_string(input_time);
time_parts=split(input_time,':');

if numel(time_parts)~=3
    error('Invalid FITS header time: %s',input_time)
end

hours_value=str2double(time_parts{1});
minutes_value=str2double(time_parts{2});
seconds_value=str2double(time_parts{3});

if ~isfinite(hours_value) || ~isfinite(minutes_value) || ~isfinite(seconds_value)
    error('Invalid FITS header time: %s',input_time)
end

if hours_value<0 || minutes_value<0 || seconds_value<0
    error('Invalid negative FITS header time: %s',input_time)
end

normalised_time_duration=seconds(hours_value*3600+minutes_value*60+seconds_value);

end

function coordinate_value=get_preferred_coordinate(headers,primary_field_name,fallback_field_name,warning_message)

%Read preferred coordinate field
coordinate_value=[];
if isfield(headers,primary_field_name) && ~isempty(headers.(primary_field_name))
    coordinate_value=headers.(primary_field_name);
elseif isfield(headers,fallback_field_name) && ~isempty(headers.(fallback_field_name))
    coordinate_value=headers.(fallback_field_name);
else
    warning(warning_message)
end

end

function output=clean_header_string(input_value)

%Clean quoted header string
if isstring(input_value), input_value=char(input_value); end
if ~ischar(input_value), output=''; return; end
output=strtrim(input_value);
output=strrep(output,'''','');
output=strrep(output,'"','');
output=strtrim(output);

end

function save_legacy_named_variable(file_name,variable_name,variable_value)

%Save MAT file with legacy variable name
save_structure=struct();
save_structure.(variable_name)=variable_value;
save(file_name,'-struct','save_structure')

end
