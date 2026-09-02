function [backfit,yy,xx,zz]=fit_background(summed_flat,allorders)
%select points for background fit

%Validate inputs
if ~isnumeric(summed_flat) || ndims(summed_flat)~=2, error('summed_flat must be a numeric 2D array.'); end
if ~isstruct(allorders) || ~isfield(allorders,'numords'), error('allorders must be a structure containing numords.'); end

%Build order-number mapping
test=95-allorders.numords;
Q=nan(allorders.numords,2);
for order_number=1:allorders.numords
    if order_number<=90-test
        mapped_order=151-order_number-test;
        Q(order_number,:)=[order_number mapped_order];
    else
        mapped_order=150-order_number-test; %due to missing order (m=60)
        Q(order_number,:)=[order_number mapped_order];
    end
end
m=Q(:,2);

%Set background sampling grid
xx=1:size(summed_flat,2);
st=160;
numords=allorders.numords;
yy=[];

%Fit inter-order background traces
trace_count=0;
for order_number=2:numords
    previous_points_field_name=['points' num2str(order_number-1)];
    current_points_field_name=['points' num2str(order_number)];
    previous_fit_field_name=['ofit' num2str(order_number-1)];
    current_fit_field_name=['ofit' num2str(order_number)];
    if ~isfield(allorders,previous_points_field_name) || ~isfield(allorders,current_points_field_name) || ~isfield(allorders,previous_fit_field_name) || ~isfield(allorders,current_fit_field_name), continue; end
    previous_points=allorders.(previous_points_field_name);
    current_points=allorders.(current_points_field_name);
    previous_fit=allorders.(previous_fit_field_name);
    current_fit=allorders.(current_fit_field_name);
    minpoint=max([min(previous_points) min(current_points)]);
    maxpoint=min([max(previous_points) max(current_points)]);
    if minpoint>=maxpoint, continue; end
    x=minpoint:st:maxpoint;
    if isempty(x), continue; end
    previous_start_index=find(previous_points==min(x),1,'first');
    current_start_index=find(current_points==min(x),1,'first');
    if isempty(previous_start_index) || isempty(current_start_index), continue; end
    previous_indices=previous_start_index:st:previous_start_index+st*(numel(x)-1);
    current_indices=current_start_index:st:current_start_index+st*(numel(x)-1);
    valid_indices=previous_indices<=numel(previous_fit) & current_indices<=numel(current_fit);
    x=x(valid_indices);
    previous_indices=previous_indices(valid_indices);
    current_indices=current_indices(valid_indices);
    if numel(x)<4, continue; end
    trace_count=trace_count+1;
    %Note from Emily- I am not sure of the function of this code- to check one day.
    if Q(order_number,2)<62
        if Q(order_number,2)>58.5 && Q(order_number,2)<59.5
            ty=diff(cat(1,previous_fit(previous_indices),current_fit(current_indices)));
            y=round(previous_fit(previous_indices)+ty/4);
        else
            y=round(mean(cat(1,previous_fit(previous_indices),current_fit(current_indices)),1));
        end
        pos=y(find_nearest(x,size(summed_flat,2)/2));
        [polynomial_coefficients,polynomial_structure,polynomial_scaling]=polyfit(x,y,3);
        current_yy=polyval(polynomial_coefficients,xx,polynomial_structure,polynomial_scaling);
        lastyy=current_yy;
        posgood=lastyy(round(numel(lastyy)/2));
        %startfrom goodind
        goodind=Q(Q(:,2)==62,1);
        if trace_count>1
            yy(trace_count,:)=yy(trace_count-1,:)+(pos-posgood);
        else
            yy(trace_count,:)=current_yy;
        end
    else
        y=round(mean(cat(1,previous_fit(previous_indices),current_fit(current_indices)),1));
        [polynomial_coefficients,polynomial_structure,polynomial_scaling]=polyfit(x,y,3);
        yy(trace_count,:)=polyval(polynomial_coefficients,xx,polynomial_structure,polynomial_scaling);
    end
end

%Stop if no background traces were fitted
if isempty(yy), error('No valid inter-order background traces were fitted.'); end

%Median-filter flat image
T=medfilt2(summed_flat,[5 5]);

%extend points to either side of end orders
if size(yy,1)>=2
    y1=yy(1,:)-diff(yy(1:2,:));
    yend=yy(end,:)+diff(yy(end-1:end,:));
else
    y1=yy(1,:);
    yend=yy(end,:);
end
yy=cat(1,y1,yy,yend);

%Sample background image at fitted traces
yy=round(yy);
zz=nan(size(yy));
for trace_index=1:size(yy,1)
    for column_index=1:size(yy,2)
        if yy(trace_index,column_index)>=1 && yy(trace_index,column_index)<=size(summed_flat,1)
            zz(trace_index,column_index)=T(yy(trace_index,column_index),xx(column_index));
        end
    end
end

%Prepare interpolation grid
xx=repmat(xx,[size(yy,1) 1]);
[Y,X]=meshgrid(1:size(summed_flat,1),1:size(summed_flat,2));

%interpolating background points
warning('off','MATLAB:griddata:DuplicateDataPoints');

%remove nans
valid_background_points=isfinite(yy) & isfinite(xx) & isfinite(zz);
if nnz(valid_background_points)<10, error('Too few valid background points for interpolation.'); end

%Interpolate background model
backfit=griddata(yy(valid_background_points),xx(valid_background_points),zz(valid_background_points),Y,X,'cubic');

%Replace missing interpolation values
B=isnan(backfit);
backfit(B)=0;
backfit=backfit';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nearest=find_nearest(vector,value)

%Find nearest index
[~,nearest]=min(abs(vector-value));

end