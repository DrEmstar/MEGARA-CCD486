function [star,background,diagnostic]=extract_registered_stellar_orders(allorders,data,backfit,flatdata,cosmic_plotting)
% Extract once to measure registration, then re-extract at the shifted trace.
star=extract_all_orders_no_background(allorders,data);
background=extract_all_orders_no_background(allorders,backfit);
star.numords=allorders.numords;
initial_plotting=0;
if isfield(allorders,'stellar_trace_registration_enabled') && ...
        ~allorders.stellar_trace_registration_enabled
    initial_plotting=cosmic_plotting;
end
star=remove_cosmics(star,initial_plotting);
[registered,diagnostic]=register_stellar_trace(allorders,star,background,flatdata);
diagnostic.fallback_orders=[];
if strcmp(diagnostic.status,'applied') && diagnostic.applied_shift_pixels~=0
    for j=1:allorders.numords
        field=['order_' num2str(j)];
        if ~isfield(star,field) || ~isfield(background,field), continue; end
        try
            points=registered.(['points' num2str(j)]);
            trace=registered.(['ofit' num2str(j)]);
            width=registered.(['width' num2str(j)]);
            factor=1.3;
            if isfield(registered,'WIDTHFACTOR') && ~isempty(registered.WIDTHFACTOR)
                factor=registered.WIDTHFACTOR;
            end
            s=extract_order_no_background(points,trace,width,data,factor);
            b=extract_order_no_background(points,trace,width,backfit,factor);
            % A detector-edge shift must not change wavelength/flat indexing.
            if ~isequal(s.xax,star.(field).xax) || ~isequal(b.xax,background.(field).xax)
                error('MEGARA:TraceGeometry','Shift changes the extracted detector columns.');
            end
            cleaned=remove_cosmics(struct('numords',1,'order_1',s),cosmic_plotting);
            star.(field)=cleaned.order_1; background.(field)=b;
        catch err
            diagnostic.fallback_orders(end+1)=j;
            warning('MEGARA:TraceOrderFallback','Order %d retained its original trace: %s',j,err.message);
        end
    end
elseif ~strcmp(diagnostic.status,'disabled') && ~strcmp(diagnostic.status,'applied')
    warning('MEGARA:TraceRegistrationFallback','Retaining monthly master trace: %s.',diagnostic.status);
end
star.cosmic_count=0;
for j=1:allorders.numords
    field=['order_' num2str(j)];
    if isfield(star,field) && isfield(star.(field),'cosmics')
        star.cosmic_count=star.cosmic_count+nnz(star.(field).cosmics);
    end
end
fprintf('Stellar trace registration: %s, shift %.3f pixels (%d orders).\n', ...
    diagnostic.status,diagnostic.applied_shift_pixels,nnz(diagnostic.accepted_orders));
end
