function [fullwave,fullint_rel]=combining_orders_2Dint(data,step,ordmin,ordmax)
%make final wavelength axis

%Check order direction
minwave1=min(data.(['wave' num2str(ordmin)]));
minwave2=min(data.(['wave' num2str(ordmax)]));

if minwave1>minwave2
    %then redo order numbering so increasing order number corresponds to increasing wavelength
    data=reverse_order_numbering(data,ordmin,ordmax);
end

%Create final wavelength grid
fullwave=rounding(min(data.(['wave' num2str(ordmin)])),step):step:rounding(max(data.(['wave' num2str(ordmax)])),step);

%Initialise output matrix
numobs=size(data.(['int' num2str(ordmin)]),1);
fullint_rel=nan(numobs,length(fullwave));

%Merge each observation
for OBS=1:numobs
    weighted_intensity_sum=zeros(1,length(fullwave));
    fullweight=zeros(1,length(fullwave));
    for order_number=ordmin:ordmax
        wave_field_name=['wave' num2str(order_number)];
        intensity_field_name=['int' num2str(order_number)];
        weights_field_name=['weights' num2str(order_number)];
        if ~isfield(data,wave_field_name) || ~isfield(data,intensity_field_name) || ~isfield(data,weights_field_name), continue; end
        wave=data.(wave_field_name);
        int=data.(intensity_field_name)(OBS,:);
        weight=data.(weights_field_name)(OBS,:);
        [wave,int,weight]=clean_order_vectors(wave,int,weight);
        if numel(wave)<2, continue; end
        weight=apply_overlap_edge_taper(data,OBS,order_number,wave,weight,ordmin,ordmax);
        order_mask=fullwave>=min(wave) & fullwave<=max(wave);
        if ~any(order_mask), continue; end
        interpolated_intensity=interp1(wave,int,fullwave(order_mask),'linear',NaN);
        interpolated_weight=interp1(wave,weight,fullwave(order_mask),'linear',NaN);
        valid_mask=isfinite(interpolated_intensity) & isfinite(interpolated_weight) & interpolated_weight>0;
        order_indices=find(order_mask);
        order_indices=order_indices(valid_mask);
        weighted_intensity_sum(order_indices)=weighted_intensity_sum(order_indices)+interpolated_intensity(valid_mask).*interpolated_weight(valid_mask);
        fullweight(order_indices)=fullweight(order_indices)+interpolated_weight(valid_mask);
    end
    valid_fullweight=fullweight>0;
    fullint_rel(OBS,valid_fullweight)=weighted_intensity_sum(valid_fullweight)./fullweight(valid_fullweight);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data=reverse_order_numbering(data,ordmin,ordmax)

%Reverse order fields
temporary_data=struct();
order_count=0;
for old_order_number=ordmin:ordmax
    order_count=order_count+1;
    new_order_number=ordmax+1-order_count;
    old_wave_field_name=['wave' num2str(old_order_number)];
    old_intensity_field_name=['int' num2str(old_order_number)];
    old_weights_field_name=['weights' num2str(old_order_number)];
    new_wave_field_name=['wave' num2str(new_order_number)];
    new_intensity_field_name=['int' num2str(new_order_number)];
    new_weights_field_name=['weights' num2str(new_order_number)];
    if isfield(data,old_wave_field_name), temporary_data.(new_wave_field_name)=data.(old_wave_field_name); end
    if isfield(data,old_intensity_field_name), temporary_data.(new_intensity_field_name)=data.(old_intensity_field_name); end
    if isfield(data,old_weights_field_name), temporary_data.(new_weights_field_name)=data.(old_weights_field_name); end
end
data=temporary_data;

end

function [wave,int,weight]=clean_order_vectors(wave,int,weight)

%Prepare vectors for interpolation
wave=wave(:)';
int=int(:)';
weight=weight(:)';
valid_wavelength=isfinite(wave);
wave=wave(valid_wavelength);
int=int(valid_wavelength);
weight=weight(valid_wavelength);
[wave,sort_index]=sort(wave);
int=int(sort_index);
weight=weight(sort_index);
[wave,unique_index]=unique(wave,'stable');
int=int(unique_index);
weight=weight(unique_index);
invalid_pixel=~isfinite(int) | ~isfinite(weight) | weight<=0;
int(invalid_pixel)=NaN;
weight(invalid_pixel)=0;

end

function weight=apply_overlap_edge_taper(data,OBS,order_number,wave,weight,ordmin,ordmax)

%Reduce the influence of order edges only where another order overlaps.
taper_length=min(50,numel(weight));
if taper_length<2, return; end

if edge_is_covered_by_neighbour(data,OBS,order_number,wave(1),ordmin,ordmax)
    weight(1:taper_length)=weight(1:taper_length).*linspace(0,1,taper_length);
end
if edge_is_covered_by_neighbour(data,OBS,order_number,wave(end),ordmin,ordmax)
    weight(end-taper_length+1:end)=weight(end-taper_length+1:end).*linspace(1,0,taper_length);
end

end

function covered=edge_is_covered_by_neighbour(data,OBS,order_number,edge_wavelength,ordmin,ordmax)

covered=false;
neighbour_orders=[order_number-1 order_number+1];
for neighbour_order=neighbour_orders
    if neighbour_order<ordmin || neighbour_order>ordmax, continue; end
    wave_field_name=['wave' num2str(neighbour_order)];
    intensity_field_name=['int' num2str(neighbour_order)];
    weights_field_name=['weights' num2str(neighbour_order)];
    if ~isfield(data,wave_field_name) || ~isfield(data,intensity_field_name) || ...
            ~isfield(data,weights_field_name), continue; end
    neighbour_wave=data.(wave_field_name);
    neighbour_intensity=data.(intensity_field_name)(OBS,:);
    neighbour_weight=data.(weights_field_name)(OBS,:);
    neighbour_wave=neighbour_wave(:)';
    neighbour_intensity=neighbour_intensity(:)';
    neighbour_weight=neighbour_weight(:)';
    if numel(neighbour_wave)~=numel(neighbour_intensity) || ...
            numel(neighbour_wave)~=numel(neighbour_weight), continue; end
    valid_neighbour=isfinite(neighbour_wave) & isfinite(neighbour_intensity) & ...
        isfinite(neighbour_weight) & neighbour_weight>0;
    if any(valid_neighbour) && edge_wavelength>=min(neighbour_wave(valid_neighbour)) && ...
            edge_wavelength<=max(neighbour_wave(valid_neighbour))
        covered=true;
        return
    end
end

end

function output=rounding(input,tonearest)

%Round to nearest step
output=round(input/tonearest)*tonearest;

end
