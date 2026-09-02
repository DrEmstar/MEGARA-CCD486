function stardata_c=remove_cosmics(stardata,plotting)

%Check number of orders
if ~isstruct(stardata) || ~isfield(stardata,'numords')
    disp('Can''t determine number of orders -> skipping')
    stardata_c=stardata;
    return
end

%Initialise counter
numords=stardata.numords;
cosmic_count=0;

%Clean each order
for order_number=1:numords
    order_field_name=['order_' num2str(order_number)];
    if ~isfield(stardata,order_field_name)
        fprintf('skipping order %.0f\n',order_number)
        continue
    end
    if ~isfield(stardata.(order_field_name),'data')
        fprintf('skipping order %.0f\n',order_number)
        continue
    end
    odata=stardata.(order_field_name).data;
    if isempty(odata) || ~isnumeric(odata)
        fprintf('skipping order %.0f\n',order_number)
        continue
    end

    %Build median-filtered replacement image
    odatab=odata;
    filtered_order_data=medfilt2(odata,[7 1]);

    %Normalise order data
    order_mean=mean(odata(:),'omitnan');
    if ~isfinite(order_mean) || order_mean==0
        fprintf('skipping order %.0f\n',order_number)
        continue
    end
    normalised_order_data=odata./order_mean;

    %Find cosmic-ray candidates
    residuals=normalised_order_data-medfilt2(normalised_order_data,[7 1]);
    cosmic_pixels=residuals>0.5;
    cosmic_count=cosmic_count+nnz(cosmic_pixels);

    %fix up odatab
    odatab(cosmic_pixels)=filtered_order_data(cosmic_pixels);

    %Store cleaned order
    stardata.(order_field_name).data=odatab;
    stardata.(order_field_name).data_pre_c=odata;
    stardata.(order_field_name).cosmics=cosmic_pixels;

    %Plot diagnostic figures
    if plotting==1
        % commands to produce test spectrum
        figure(100)
        mesh(odatab)
        title(['cosmic ray fixed, order' num2str(order_number)])
        figure(101)
        mesh(odata)
        title(['before cosmic ray fix, order' num2str(order_number)])
        pause
    end
end

%Store output
stardata.cosmic_count=cosmic_count;
stardata_c=stardata;

end