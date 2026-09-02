function [data2,datacf]=harmonise_continuum_orders( ...
    data2,datacf,firstorder,lastorder,step)
%Match adjacent order levels using small regularised scale corrections.

maximum_adjustment=0.05;
regularisation_strength=0.10;

if isfield(datacf,'automatic_order_scales')
    orders=datacf.automatic_order_scales.orders;
    scale_factors=datacf.automatic_order_scales.scale_factors;
else
    [orders,scale_factors,overlap_offsets]=calculate_order_scales( ...
        data2,firstorder,lastorder,step,regularisation_strength, ...
        maximum_adjustment);
    datacf.automatic_order_scales.orders=orders;
    datacf.automatic_order_scales.scale_factors=scale_factors;
    datacf.automatic_order_scales.overlap_offsets=overlap_offsets;
    datacf.automatic_order_scales.maximum_adjustment=maximum_adjustment;
    datacf.automatic_order_scales.regularisation_strength= ...
        regularisation_strength;
end

for order_index=1:numel(orders)
    intensity_field=['int' num2str(orders(order_index))];
    if ~isfield(data2,intensity_field), continue; end
    data2.(intensity_field)=data2.(intensity_field).* ...
        scale_factors(order_index);
end

end

function [orders,scale_factors,overlap_offsets]=calculate_order_scales( ...
    data2,firstorder,lastorder,step,regularisation_strength, ...
    maximum_adjustment)

orders=[];
for order_number=firstorder:lastorder
    if isfield(data2,['wave' num2str(order_number)]) && ...
            isfield(data2,['int' num2str(order_number)])
        orders(end+1)=order_number; %#ok<AGROW>
    end
end

number_orders=numel(orders);
overlap_matrix=[];
overlap_target=[];
overlap_rows=cell(0,4);

for order_index=1:number_orders-1
    first_order=orders(order_index);
    second_order=orders(order_index+1);
    if second_order~=first_order+1, continue; end

    first_wave=data2.(['wave' num2str(first_order)]);
    second_wave=data2.(['wave' num2str(second_order)]);
    first_intensity=median(data2.(['int' num2str(first_order)]), ...
        1,'omitnan');
    second_intensity=median(data2.(['int' num2str(second_order)]), ...
        1,'omitnan');
    overlap_start=max(min(first_wave),min(second_wave));
    overlap_stop=min(max(first_wave),max(second_wave));
    if overlap_start>=overlap_stop, continue; end

    overlap_wave=overlap_start:step:overlap_stop;
    first=interp1(first_wave,first_intensity,overlap_wave,'linear',NaN);
    second=interp1(second_wave,second_intensity,overlap_wave,'linear',NaN);
    valid=isfinite(first) & isfinite(second) & first>0 & second>0;
    if nnz(valid)<20, continue; end

    log_ratio=median(log(first(valid)./second(valid)),'omitnan');
    equation=zeros(1,number_orders);
    equation(order_index)=1;
    equation(order_index+1)=-1;
    overlap_matrix(end+1,:)=equation; %#ok<AGROW>
    overlap_target(end+1,1)=-log_ratio; %#ok<AGROW>
    overlap_rows(end+1,:)={first_order,second_order,nnz(valid), ...
        exp(log_ratio)-1}; %#ok<AGROW>
end

regularisation=sqrt(regularisation_strength)*eye(number_orders);
solution=[overlap_matrix;regularisation]\ ...
    [overlap_target;zeros(number_orders,1)];
minimum_log_scale=log(1-maximum_adjustment);
maximum_log_scale=log(1+maximum_adjustment);
solution=min(max(solution,minimum_log_scale),maximum_log_scale);
scale_factors=exp(solution).';
overlap_offsets=cell2table(overlap_rows,'VariableNames', ...
    {'FirstOrder','SecondOrder','NumberPixels','FractionalOffset'});

fprintf(['Adjusted the relative continuum levels of %.0f orders to improve ' ...
    'agreement in their overlap regions\n'],number_orders)

end
