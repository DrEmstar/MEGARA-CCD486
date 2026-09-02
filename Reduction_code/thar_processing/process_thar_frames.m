function [jd,m,fwave,numbad_thar_flag,info_rebuild,kept_line_count,instrument_lines]=process_thar_frames(tharname,info,allorders,Bias,flatdata,blue_data_chop_value,move_rejected_frame)

%Initialise outputs
jd=NaN;
m=[];
fwave=[];
numbad_thar_flag=0;
instrument_lines=table();
info_rebuild=0;
kept_line_count=NaN;
if nargin<7, move_rejected_frame=true; end

%Resolve ThAr path and display name
[thar_folder,thar_base,thar_ext]=fileparts(tharname);
if isempty(thar_ext)
    thar_ext='.fit';
end
if isempty(thar_folder)
    thar_file_path=fullfile(pwd,[thar_base thar_ext]);
    thar_folder=pwd;
else
    thar_file_path=tharname;
end
thar_display_name=[thar_base thar_ext];

%ThAr image processing (1 image at a time only)
backfit.back_polyfit=0; %ThAr images don't have background measurements
backfit.xscale=1;backfit.xcen=0;
backfit.yscale=1;backfit.ycen=0;
backfit.zscale=1;backfit.zcen=0;

%Read observation metadata
sname=thar_base;
if ~isfield(info,sname), error(['info does not contain metadata for ' sname '.']); end
if ~isfield(info.(sname),'MJD_OBS'), error(['info.' sname ' does not contain MJD_OBS.']); end
jd=info.(sname).MJD_OBS;

%Read ThAr image
data=raw_hercules_fitsread(thar_file_path);

%assuming no Bias image
if ~isempty(Bias)
    if isequal(size(Bias),size(data))
        data=data-Bias;
    else
        warning(['Bias size does not match ' thar_display_name '. Bias subtraction skipped.'])
    end
end

%Chop blue detector rows
data=blue_data_chop(data,blue_data_chop_value);

%Extract ThAr orders
thardata=extract_all_orders_no_background(allorders,data);

%Apply flat-profile weighted extraction
thardata=apply_flat_profile_weights_to_thorium(thardata,flatdata);

%Handling rollover of 4 digit filenames from J5XXXX to J6XXXX
thardate=str2double(sname(2:5));
if ~isfinite(thardate), error(['Could not parse ThAr date code from filename: ' thar_display_name]); end

%Load ThAr wavelength reference data
[order,air,ul,uc,ur]=load_thorium_reference_data(thardate,jd);

%apply simple program to change to absolute order numbers
[thardata_2,m]=redo_order_numbering_simple(thardata,allorders);

%finding and fitting ThAr lines in extracted spectrum
%make plot_q=1 if you want plots of every th line as it is processed to
%display on screen pausing for the enter button to continue
plot_q=0;

%Find ThAr lines quietly
captured_find_output=evalc('[thfitinfo,numbad_thar_flag,move_file_flag]=find_th_lines(thardata_2,order,ul,uc,ur,plot_q);');
kept_line_count=parse_kept_thar_line_count(captured_find_output);

%Handle manually rejected ThAr frame
if move_file_flag==1

    if move_rejected_frame
        %Create the manual_removed folder in this raw folder
        targetFolder=fullfile(thar_folder,'manual_removed');
        if ~exist(targetFolder,'dir'), mkdir(targetFolder); end

        %Move the rejected ThAr file using an absolute path
        [move_status,move_message]=movefile(thar_file_path,targetFolder,'f');
        if ~move_status
            warning(['Could not move manually rejected ThAr file ' thar_display_name ': ' move_message])
        end
    end

    info_rebuild=1;
    fwave=[];

else

    %wavelength calibration
    captured_wavelength_output=evalc('wavfit=wavelength_calibration(thfitinfo,order,air);');

    if isnan(kept_line_count)
        kept_line_count=parse_kept_thar_line_count(captured_wavelength_output);
    end

    %make wavelength axes for each ThAr for data:
    fwave=make_wave_vector(thardata_2,wavfit,m);
    instrument_lines=measure_instrument_lines(thfitinfo,order,thardata_2,fwave);
    info_rebuild=0;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kept_line_count=parse_kept_thar_line_count(captured_output)

%Default if line count cannot be parsed
kept_line_count=NaN;

%Parse common retained-line output formats
tokens=regexp(captured_output,'(\d+)\s+lines\s+chosen','tokens');

if isempty(tokens)
    tokens=regexp(captured_output,'(\d+)\s+lines\s+kept','tokens');
end

if isempty(tokens)
    tokens=regexp(captured_output,'(\d+)\s+lines\s+used','tokens');
end

if ~isempty(tokens)
    kept_line_count=str2double(tokens{end}{1});
end

end

function thardata=apply_flat_profile_weights_to_thorium(thardata,flatdata)

%Check flat order count
if ~isfield(flatdata,'numords'), error('flatdata does not contain numords.'); end

%Apply extraction profile to each order
for order_number=1:flatdata.numords
    order_field_name=['order_' num2str(order_number)];
    if ~isfield(flatdata,order_field_name), warning(['Skipping ThAr order ' num2str(order_number) ' because flatdata.' order_field_name ' is missing.']); continue; end
    if ~isfield(thardata,order_field_name), warning(['Skipping ThAr order ' num2str(order_number) ' because thardata.' order_field_name ' is missing.']); continue; end
    if ~isfield(flatdata.(order_field_name),'extraction_inds'), error(['flatdata.' order_field_name ' does not contain extraction_inds.']); end
    if ~isfield(thardata.(order_field_name),'data'), error(['thardata.' order_field_name ' does not contain data.']); end
    ind=flatdata.(order_field_name).extraction_inds;
    if isfield(flatdata.(order_field_name),'profile_weights')
        profile=flatdata.(order_field_name).profile_weights;
    else
        profile=zeros(size(ind));
        profile(ind)=1;
        profile_sum=sum(profile);
        if profile_sum==0 || ~isfinite(profile_sum), profile=ones(size(ind)); profile_sum=sum(profile); end
        profile=profile./profile_sum;
    end
    odata=thardata.(order_field_name).data;
    if numel(ind)~=size(odata,2), error(['Extraction index length does not match ThAr data width for ' order_field_name '.']); end
    if numel(profile)~=size(odata,2), error(['Profile weight length does not match ThAr data width for ' order_field_name '.']); end
    thardata.(order_field_name).summed_data=sum(odata(:,ind).*profile(ind),2);
end

end

function [order,air,ul,uc,ur]=load_thorium_reference_data(thardate,jd)

%Select ThAr reference file
thorium_reference_file_name=select_thorium_reference_file(thardate,jd);

%Load MAT reference file if available
thorium_reference_file_path=find_calibration_file(thorium_reference_file_name);
if ~isempty(thorium_reference_file_path)
    loaded_reference=load(thorium_reference_file_path);
    [order,air,ul,uc,ur,reference_loaded]=extract_thorium_reference_variables(loaded_reference);
    if reference_loaded
        return
    end
end

%Fallback to thar.dat
[order,air,ul,uc,ur]=read_thar_dat_file();

end

function calibration_file_path=find_calibration_file(calibration_file_name)

%Find calibration file in current folder or MATLAB path
calibration_file_path='';
if isempty(calibration_file_name), return; end
if isfile(calibration_file_name)
    calibration_file_path=calibration_file_name;
    return
end
path_match=which(calibration_file_name);
if ~isempty(path_match)
    calibration_file_path=path_match;
end

end

function thorium_reference_file_name=select_thorium_reference_file(thardate,jd)

%Select historical ThAr reference database
thorium_reference_file_name='';
if jd<60000
    %Dynamic thardat selection tested for 2007-2026 data
    if thardate<=4463
        thorium_reference_file_name='thardat_g.mat';
    elseif thardate>4463 && thardate<=4781
        thorium_reference_file_name='thardat_f.mat';
    elseif thardate>4781 && thardate<=5454
        thorium_reference_file_name='thardat_d.mat';
    elseif thardate>5454 && thardate<=5702
        thorium_reference_file_name='thardat_e.mat';
    elseif thardate>5702 && thardate<=6337
        thorium_reference_file_name='thardat_c.mat';
    elseif thardate>6337 && thardate<6791
        thorium_reference_file_name='thardat_b.mat';
    elseif thardate>6791 && thardate<7083
        thorium_reference_file_name='thardat_a.mat';
    elseif thardate>=7083 && thardate<7178
        thorium_reference_file_name='thardat_i.mat';
    elseif thardate>=7178 && thardate<7821
        thorium_reference_file_name='thardat_a.mat';
    elseif thardate>=7821 && thardate<7968
        thorium_reference_file_name='thardat_h.mat';
    elseif thardate>=7968 && thardate<8120
        thorium_reference_file_name='thardat_a.mat';
    elseif thardate>=8120 && thardate<8144
        thorium_reference_file_name='thardat_j.mat';
    elseif thardate>=8144 && thardate<8146
        thorium_reference_file_name='thardat_k.mat';
    elseif thardate>=8146 && thardate<8303
        thorium_reference_file_name='thardat_j.mat';
    elseif thardate>=8303 && thardate<9461
        thorium_reference_file_name='thardat_l.mat';
    elseif thardate>=9461 && thardate<9650
        thorium_reference_file_name='thardat_mm.mat';
    elseif thardate>=9650
        thorium_reference_file_name='thardat_n.mat';
    end
elseif thardate<889
    thorium_reference_file_name='thardat_o.mat';
else
    thorium_reference_file_name='thardat_p.mat';
end

end

function [order,air,ul,uc,ur,reference_loaded]=extract_thorium_reference_variables(loaded_reference)

%Initialise outputs
order=[];
air=[];
ul=[];
uc=[];
ur=[];
reference_loaded=false;

%Extract direct variables
if isfield(loaded_reference,'order') && isfield(loaded_reference,'air') && isfield(loaded_reference,'ul') && isfield(loaded_reference,'uc') && isfield(loaded_reference,'ur')
    order=loaded_reference.order;
    air=loaded_reference.air;
    ul=loaded_reference.ul;
    uc=loaded_reference.uc;
    ur=loaded_reference.ur;
    reference_loaded=true;
    return
end

%Extract variables from first contained structure if needed
loaded_field_names=fieldnames(loaded_reference);
for field_index=1:numel(loaded_field_names)
    candidate=loaded_reference.(loaded_field_names{field_index});
    if isstruct(candidate) && isfield(candidate,'order') && isfield(candidate,'air') && isfield(candidate,'ul') && isfield(candidate,'uc') && isfield(candidate,'ur')
        order=candidate.order;
        air=candidate.air;
        ul=candidate.ul;
        uc=candidate.uc;
        ur=candidate.ur;
        reference_loaded=true;
        return
    end
end

end

function [order,air,ul,uc,ur]=read_thar_dat_file()

%Read fallback ThAr text file
thar_dat_file=find_calibration_file('thar.dat');
if isempty(thar_dat_file)
    error('there is no thardat_*.mat or thar.dat file accessible; add the folder containing the ThAr calibration files to the MATLAB path, or copy the required files into the current reduction folder.')
end
file_identifier=fopen(thar_dat_file,'r');
if file_identifier==-1, error(['Could not open thar.dat: ' thar_dat_file]); end
cleanup_file=onCleanup(@() fclose(file_identifier));
C=textscan(file_identifier,'%.0f %.4f %.4f %.4f %s %.3f 2017%.3f %.3f %.3f %.3f %.3f');
order=C{1};
air=C{2};
vac=C{3};
wavnum=C{4};
species=C{5};
ul=C{6};
uc=C{7};
ur=C{8};
vl=C{9};
vc=C{10};
vr=C{11};

% standard corrections to the text file
uc=uc/0.015+36;vc=vc/0.015-200;
ul=ul/0.015+36;vl=vl/0.015-200;
ur=ur/0.015+36;vr=vr/0.015-200;

end

function instrument_lines=measure_instrument_lines(thfitinfo,line_orders,thardata,fwave)

%Convert accepted ThAr Gaussian widths from pixels to velocity
c=299792.458;
number_lines=size(thfitinfo,1);
wavelength=nan(number_lines,1);
fwhm_pixels=thfitinfo(:,5);
fwhm_velocity=nan(number_lines,1);
order_number=nan(number_lines,1);
fit_resnorm=thfitinfo(:,2);

for fit_index=1:number_lines
    reference_index=round(thfitinfo(fit_index,1));
    if reference_index<1 || reference_index>numel(line_orders), continue; end
    order_number(fit_index)=line_orders(reference_index);
    order_field_name=['order_' num2str(order_number(fit_index))];
    if ~isfield(thardata,order_field_name) || ~isfield(fwave,order_field_name), continue; end
    pixel=thardata.(order_field_name).xax(:);
    order_wavelength=fwave.(order_field_name);
    order_wavelength=order_wavelength(:);
    if numel(pixel)~=numel(order_wavelength) || numel(pixel)<3, continue; end
    line_position=thfitinfo(fit_index,4);
    wavelength(fit_index)=interp1(pixel,order_wavelength,line_position,'linear',NaN);
    wavelength_gradient=gradient(order_wavelength,pixel);
    local_dispersion=interp1(pixel,wavelength_gradient,line_position,'linear',NaN);
    if isfinite(wavelength(fit_index)) && wavelength(fit_index)>0 && ...
            isfinite(local_dispersion)
        fwhm_velocity(fit_index)=c*fwhm_pixels(fit_index)* ...
            abs(local_dispersion)/wavelength(fit_index);
    end
end

valid=isfinite(wavelength) & isfinite(fwhm_velocity) & ...
    fwhm_velocity>0 & fwhm_velocity<100;
instrument_lines=table(wavelength(valid),order_number(valid), ...
    fwhm_pixels(valid),fwhm_velocity(valid),fit_resnorm(valid), ...
    'VariableNames',{'Wavelength','Order','FwhmPixels','FwhmVelocity', ...
    'FitResnorm'});

end
