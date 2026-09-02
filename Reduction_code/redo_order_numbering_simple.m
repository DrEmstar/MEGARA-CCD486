function [Data,m]=redo_order_numbering_simple(data,allorders)

%Initialise output
Data=struct();
m=[];

%Check order information
if ~isstruct(allorders) || ~isfield(allorders,'numords'), error('allorders must contain numords.'); end
number_of_orders=allorders.numords;
if ~(isnumeric(number_of_orders) && isscalar(number_of_orders) && isfinite(number_of_orders))
    error('allorders.numords must be a finite numeric scalar. Current class is %s.',class(number_of_orders))
end
number_of_orders=double(number_of_orders);
if number_of_orders<2, error('At least two traced orders are required for order renumbering. allorders.numords=%.0f.',number_of_orders); end

%find the mean y-pixel position of each of the orders
meanordypos=nan(number_of_orders,1);
for order_number=1:number_of_orders
    order_fit_field_name=['ofit' num2str(order_number)];
    if ~isfield(allorders,order_fit_field_name), error(['allorders is missing ' order_fit_field_name '.']); end
    meanordypos(order_number)=mean(allorders.(order_fit_field_name),'omitnan');
end

%find where the found orders skips one, that skipped order is order 60
if number_of_orders<70, error('Order renumbering assumes at least 70 traced orders. Current allorders.numords=%.0f.',number_of_orders); end
order_spacing=diff(meanordypos(70:end));
if isempty(order_spacing) || all(~isfinite(order_spacing)), error('Could not identify missing-order gap for order renumbering.'); end
[~,gap_index]=max(order_spacing);
ind=gap_index+70; %the order (in initial numbering) that indicates order 59

%Build absolute order-number mapping
m=nan(1,number_of_orders);
m(1:ind-1)=61+(ind-1)-1:-1:61;
m(ind:number_of_orders)=59:-1:59-(number_of_orders-ind);

%Copy orders into absolute order-number fields
for order_number=1:number_of_orders
    input_order_field_name=['order_' num2str(order_number)];
    output_order_field_name=['order_' num2str(m(order_number))];
    if ~isfield(data,input_order_field_name), warning(['Skipping missing input order: data.' input_order_field_name]); continue; end
    Data.(output_order_field_name)=data.(input_order_field_name);
end

end