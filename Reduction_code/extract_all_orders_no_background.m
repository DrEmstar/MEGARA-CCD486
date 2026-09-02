function extracted=extract_all_orders_no_background(allorders,img)

%Initialise output
extracted=struct();

%get order numbers
if ~isstruct(allorders) || ~isfield(allorders,'numords')
    disp('order numbers not defined in allorders.numords ... exiting')
    extracted=[];
    return
end
orders=1:allorders.numords;

%determine WIDTHFACTOR
if isfield(allorders,'WIDTHFACTOR') && ~isempty(allorders.WIDTHFACTOR)
    widthfactor=allorders.WIDTHFACTOR;
else
    widthfactor=1.3;
end

%sequentially extract each order
for order_number=orders
    %get the things needed for each order - if they're not all there then skip the order
    points_field_name=['points' num2str(order_number)];
    order_fit_field_name=['ofit' num2str(order_number)];
    width_field_name=['width' num2str(order_number)];
    if ~isfield(allorders,points_field_name) || ~isfield(allorders,order_fit_field_name) || ~isfield(allorders,width_field_name)
        fprintf('\nSkipping order %d as fitting information is not complete \n\n',order_number)
        continue
    end
    points=allorders.(points_field_name);
    ofit=allorders.(order_fit_field_name);
    width=allorders.(width_field_name);
    try
        %extract the order
        extracted_order=extract_order_no_background(points,ofit,width,img,widthfactor);
        %add extracted order to the 'extracted' output structure
        extracted.(['order_' num2str(order_number)])=extracted_order;
    catch extraction_error
        fprintf('order %d has a problem and will not be extracted\n',order_number)
        warning(extraction_error.message)
    end
end

end