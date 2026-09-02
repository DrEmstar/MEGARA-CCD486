function [allorders,flatdata,backdata,summed_flat,backfit,save_flag]=flat_processing(info,Bias,blue_data_chop_value,backdata_medfilt,plotting)

%Flat image processing
save_flag=1;

%Build cached flat filename
flatdataname=['flat_info_blue_' num2str(blue_data_chop_value) '.mat'];
legacy_flatdataname=['flat_info_blue_' num2str(blue_data_chop_value) ...
    '_Megara_1.6.mat'];
loaded_flat_information=[];
number_white_lamps=count_white_lamp_frames(info);

%A compatible MEGARA 2.1 cache always takes precedence, even when raw flats
%are still present. A legacy cache is rebuilt from raw flats when possible
%and is upgraded only when the observing run contains no raw flats.
if isfile(flatdataname)
    loaded_flat_information=load(flatdataname);
    validate_cached_flat_information(loaded_flat_information,flatdataname)
    if isfield(loaded_flat_information.allorders,'WIDTHFACTOR')
        allorders=loaded_flat_information.allorders;
        flatdata=loaded_flat_information.flatdata;
        backdata=loaded_flat_information.backdata;
        summed_flat=loaded_flat_information.summed_flat;
        backfit=loaded_flat_information.backfit;
        disp('Compatible MEGARA 2.1 flat cache loaded; rebuild not required.')
        save_flag=0;
        return
    end

    fprintf('Cached flat file %s is from MEGARA version 1.6.\n',flatdataname)
    if number_white_lamps>0
        backup_flatdataname=create_legacy_flat_cache_name( ...
            blue_data_chop_value);
        fprintf(['Raw White L flats are available; preserving the legacy ' ...
            'cache as %s and rebuilding for MEGARA 2.1.\n'], ...
            backup_flatdataname)
        movefile(flatdataname,backup_flatdataname,'f');
        loaded_flat_information=[];
    end
end

if isempty(loaded_flat_information) && number_white_lamps>0
    fprintf(['No compatible MEGARA 2.1 flat cache found. Rebuilding from ' ...
        '%d raw White L flat frames.\n'],number_white_lamps)
elseif isempty(loaded_flat_information) && isfile(legacy_flatdataname)
    loaded_flat_information=load(legacy_flatdataname);
    validate_cached_flat_information(loaded_flat_information,legacy_flatdataname)
    fprintf(['No raw White L flats or MEGARA 2.1 cache found; upgrading ' ...
        'legacy flat cache %s.\n'],legacy_flatdataname)
end

%Rebuild the MEGARA 2.1 extraction products from a legacy summed flat and
%background model only when no raw White L frames are available.
if ~isempty(loaded_flat_information)
    [allorders,flatdata,backdata,summed_flat,backfit]= ...
        upgrade_legacy_flat_information(loaded_flat_information, ...
        backdata_medfilt,plotting);
    disp('Legacy flat cache upgraded; the MEGARA 2.1 cache will be saved.')
    return
end

%Merge flat frames
summed_flat=merge_flats(info,Bias);

%Chop blue end of detector
summed_flat=blue_data_chop(summed_flat,blue_data_chop_value);

%Trace orders
minpixseparation=10;
[allorders,summed_flat]=trace_orders(summed_flat,minpixseparation);

%Set extraction aperture width factor for flat, background, ThAr, and stellar extraction
allorders.WIDTHFACTOR=3;

%Fit flat background
disp('fitting the background ...')
backfit=fit_background(summed_flat,allorders);

%Subtract flat background
summed_flat_b=summed_flat-backfit;

%Extract flat data
disp('extracting the data ...')
flatdata=extract_all_orders_no_background(allorders,summed_flat_b);

%Extract background data
disp('extracting the background data ...')
backdata=extract_all_orders_no_background(allorders,backfit);

%cosmic ray filtering, just in case
flatdata.numords=allorders.numords;
flatdata=remove_cosmics(flatdata,plotting);

[flatdata,backdata]=finalise_flat_field_data( ...
    flatdata,backdata,allorders,backdata_medfilt);

%Report completion
disp('Flat-field images processing complete')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_cached_flat_information(loaded_flat_information,cache_name)

required_variables={'allorders','flatdata','backdata','summed_flat','backfit'};
for required_variable_index=1:numel(required_variables)
    required_variable=required_variables{required_variable_index};
    if ~isfield(loaded_flat_information,required_variable)
        error('Cached flat file %s is missing variable %s.', ...
            cache_name,required_variable)
    end
end

end

function number_white_lamps=count_white_lamp_frames(info)

number_white_lamps=0;
if ~isstruct(info) || ~isfield(info,'list')
    error('info must be a structure containing info.list.')
end
for file_index=1:numel(info.list)
    file_identifier=info.list{file_index};
    if isfield(info,file_identifier) && ...
            isfield(info.(file_identifier),'HERCEXPT') && ...
            strcmp(info.(file_identifier).HERCEXPT,'White L')
        number_white_lamps=number_white_lamps+1;
    end
end

end

function backup_flatdataname=create_legacy_flat_cache_name(blue_data_chop_value)

backup_flatdataname=['flat_info_blue_' num2str(blue_data_chop_value) ...
    '_Megara_1.6.mat'];
if isfile(backup_flatdataname)
    timestamp=datestr(now,'yyyymmdd_HHMMSS');
    backup_flatdataname=['flat_info_blue_' num2str(blue_data_chop_value) ...
        '_Megara_1.6_' timestamp '.mat'];
end

end

function [allorders,flatdata,backdata,summed_flat,backfit]= ...
    upgrade_legacy_flat_information(loaded_flat_information, ...
    backdata_medfilt,plotting)

allorders=loaded_flat_information.allorders;
summed_flat=loaded_flat_information.summed_flat;
backfit=loaded_flat_information.backfit;

if isempty(summed_flat) || isempty(backfit) || ...
        ~isequal(size(summed_flat),size(backfit))
    error(['Legacy flat cache does not contain compatible summed-flat ' ...
        'and background images.'])
end

%Re-extract with the MEGARA 2.1 aperture and construct its spatial-profile
%weights from the cached calibration images.
allorders.WIDTHFACTOR=3;
summed_flat_b=summed_flat-backfit;
disp('extracting the cached flat data ...')
flatdata=extract_all_orders_no_background(allorders,summed_flat_b);
disp('extracting the cached flat background ...')
backdata=extract_all_orders_no_background(allorders,backfit);
flatdata.numords=allorders.numords;
flatdata=remove_cosmics(flatdata,plotting);
[flatdata,backdata]=finalise_flat_field_data( ...
    flatdata,backdata,allorders,backdata_medfilt);

end

function [flatdata,backdata]=finalise_flat_field_data( ...
    flatdata,backdata,allorders,backdata_medfilt)

disp('finalising flat-field data')
for order_number=1:allorders.numords
    order_field_name=['order_' num2str(order_number)];
    if ~isfield(flatdata,order_field_name), error(['flatdata is missing ' order_field_name '.']); end
    if ~isfield(backdata,order_field_name), error(['backdata is missing ' order_field_name '.']); end
    if ~isfield(flatdata.(order_field_name),'data'), error(['flatdata.' order_field_name ' is missing data.']); end
    if ~isfield(backdata.(order_field_name),'data'), error(['backdata.' order_field_name ' is missing data.']); end
    odata=flatdata.(order_field_name).data;

    %Profile-weighted extraction weights from the flat spatial profile
    profile=median(odata,1,'omitnan');
    profile(~isfinite(profile))=0;
    profile(profile<0)=0;
    profile_sum=sum(profile);
    if profile_sum==0 || ~isfinite(profile_sum)
        profile=ones(size(profile));
        profile_sum=sum(profile);
    end
    profile=profile./profile_sum;
    ind=profile>max(profile)/100;
    if ~any(ind), ind=true(size(profile)); end

    flatdata.(order_field_name).extraction_inds=ind;
    flatdata.(order_field_name).profile_weights=profile;
    flatdata.(order_field_name).sumdata= ...
        sum(odata(:,ind).*profile(ind),2)';

    bdata=backdata.(order_field_name).data;
    backdata.(order_field_name).summed_data= ...
        sum(bdata(:,ind).*profile(ind),2);
    backdata.(order_field_name).summed_smoothed_data= ...
        mean_smoothing(backdata.(order_field_name).summed_data, ...
        backdata_medfilt);
    flatdata.(order_field_name).sumdata_b= ...
        flatdata.(order_field_name).sumdata(:)- ...
        backdata.(order_field_name).summed_smoothed_data(:);
end


end
