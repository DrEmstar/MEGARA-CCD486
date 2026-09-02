function [allorders,summed_flat]=trace_orders(summed_flat,minpixseparation)

%Validate inputs
if ~isnumeric(summed_flat) || ndims(summed_flat)~=2, error('summed_flat must be a numeric 2D array.'); end
if ~isnumeric(minpixseparation) || ~isscalar(minpixseparation) || ~isfinite(minpixseparation) || minpixseparation<1, error('minpixseparation must be a finite scalar value >= 1.'); end

%Initialise output
allorders=struct();
allorders.numords=0;

%Build central spatial profile
image_size=size(summed_flat);
spatial_pixels=1:image_size(1);
centre_column=round(image_size(2)/2);
column_window=max(1,centre_column-5):min(image_size(2),centre_column+5);
tracedata=cat(2,spatial_pixels(:),mean(summed_flat(:,column_window),2,'omitnan'));
tracedata(:,2)=tracedata(:,2)-min(tracedata(:,2));

%Remove smooth background profile
smooth_background=running_min(tracedata(:,2),10);
smooth_background=mean_smoothing(medfilt1(smooth_background,21),10);
tracedata(:,2)=tracedata(:,2)-smooth_background;

%Normalise order profile
normalisation_profile=mean_smoothing(running_max(tracedata(:,2),30),30);
normalisation_profile(normalisation_profile==0 | ~isfinite(normalisation_profile))=NaN;
tracedata(:,2)=tracedata(:,2)./normalisation_profile;

%Find candidate order blocks
orders=find_order_block_centres(tracedata,minpixseparation);
if isempty(orders), warning('No candidate order centres found in trace_orders.'); return; end

disp('getting starting centre orders ...')

%fit 1-D Gaussian at 20 pixels around original trace
width=0.4;
fitsize=15;
centre=zeros(numel(orders),1);
widths=centre;
sse=centre;
for order_index=1:numel(orders)
    fit_start=max(1,round(orders(order_index)-fitsize));
    fit_stop=min(size(tracedata,1),round(orders(order_index)+fitsize));
    fitdata=tracedata(fit_start:fit_stop,:);
    if size(fitdata,1)<5
        centre(order_index)=NaN;
        widths(order_index)=NaN;
        sse(order_index)=Inf;
        continue
    end
    [centre(order_index),widths(order_index),sse(order_index),fit_values]=fit_lin_gauss(fitdata(:,1),fitdata(:,2));
end

%Reject poor starting fits
good_starting_fit=sse<0.2 & abs(widths)<10 & abs(widths)>2.5 & isfinite(centre) & isfinite(widths);
centre=centre(good_starting_fit);
widths=widths(good_starting_fit);
sse=sse(good_starting_fit);

%now we have starting centres for orders
%run an order tracing program that finds the order starting from the centre

disp('tracing each order ...')

%Trace each order
order_count=0;
for centre_index=1:numel(centre)
    [opos,ocentres,owidths]=trace_single_order_from_centre(summed_flat,centre(centre_index),widths(centre_index));
    bad_widths=owidths<2.5 | owidths>9 | ~isfinite(owidths) | ~isfinite(ocentres) | ~isfinite(opos);
    opos(bad_widths)=[];
    ocentres(bad_widths)=[];
    owidths(bad_widths)=[];
    if numel(opos)>10
        order_count=order_count+1;
        points=min(opos):max(opos);
        [polynomial_coefficients,polynomial_structure,polynomial_scaling]=polyfit(opos,ocentres,5);
        [ofit,delta]=polyval(polynomial_coefficients,points,polynomial_structure,polynomial_scaling);
        allorders.(['width' num2str(order_count)])=median(owidths);
        allorders.(['ofit' num2str(order_count)])=ofit;
        allorders.(['points' num2str(order_count)])=points;
    else
    end
end

%Store number of traced orders
allorders.numords=order_count;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function orders=find_order_block_centres(tracedata,minpixseparation)

%Find thresholded order pixels
threshold_pixels=find(tracedata(:,2)>0.25);
if isempty(threshold_pixels)
    orders=[];
    return
end

%Find separated threshold blocks
block_breaks=find(diff(threshold_pixels)>minpixseparation);
block_starts=[threshold_pixels(1); threshold_pixels(block_breaks+1)];
block_stops=[threshold_pixels(block_breaks); threshold_pixels(end)];

%Discard first and last candidate blocks as in original logic
if numel(block_starts)<=2
    orders=[];
    return
end
block_starts=block_starts(2:end-1);
block_stops=block_stops(2:end-1);

%Calculate block centres
orders=mean(cat(2,block_starts(:),block_stops(:)),2);

end