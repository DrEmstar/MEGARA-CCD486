function [registered,diagnostic]=register_stellar_trace(allorders,star,background,flatdata)
%REGISTER_STELLAR_TRACE Measure a robust science-minus-monthly-flat offset.
% Input profiles must share detector columns and trace-relative apertures.
% Only the returned science trace is changed; the master flat is immutable.
registered=allorders;
N=allorders.numords;
diagnostic=struct('method','monthly-flat-centroid-v1','enabled',true, ...
    'status','insufficient_orders','measured_shift_pixels',NaN, ...
    'applied_shift_pixels',0,'order_shift_pixels',nan(1,N), ...
    'valid_columns',zeros(1,N),'accepted_orders',false(1,N), ...
    'order_scatter_pixels',NaN);
if isfield(allorders,'stellar_trace_registration_enabled') && ...
        ~allorders.stellar_trace_registration_enabled
    diagnostic.enabled=false; diagnostic.status='disabled'; return
end
for j=1:N
    field=['order_' num2str(j)];
    if ~isfield(star,field) || ~isfield(background,field) || ~isfield(flatdata,field), continue; end
    s=star.(field); b=background.(field); f=flatdata.(field);
    if ~isfield(f,'extraction_inds') || ~isfield(f,'xax') || ...
            ~isequal(s.xax(:),f.xax(:)) || ~isequal(s.xax(:),b.xax(:)) || ...
            ~isequal(size(s.data),size(f.data),size(b.data))
        continue
    end
    ind=f.extraction_inds;
    if ~islogical(ind) || numel(ind)~=size(s.data,2) || nnz(ind)<3, continue; end
    offsets=s.yax-s.ypositions(:);
    source=s.data-b.data; master=f.data;
    good=all(isfinite(source(:,ind)),2)&all(isfinite(master(:,ind)),2);
    source=max(source(:,ind),0); master=max(master(:,ind),0);
    sf=sum(source,2); ff=sum(master,2);
    if ~any(good&sf>0&ff>0), continue; end
    sm=median(sf(good&sf>0)); fm=median(ff(good&ff>0));
    good=good&sf>.3*sm&ff>.3*fm;
    delta=sum(source.*offsets(:,ind),2)./sf- ...
        sum(master.*offsets(:,ind),2)./ff;
    good=good&isfinite(delta);
    if nnz(good)<100, continue; end
    centre=median(delta(good)); scatter=1.4826*median(abs(delta(good)-centre));
    good=good&abs(delta-centre)<=max(.25,4*scatter);
    diagnostic.valid_columns(j)=nnz(good);
    if nnz(good)<100 || scatter>1, continue; end
    diagnostic.order_shift_pixels(j)=median(delta(good));
end
valid=isfinite(diagnostic.order_shift_pixels);
if nnz(valid)<3, return; end
shift=median(diagnostic.order_shift_pixels(valid));
scatter=1.4826*median(abs(diagnostic.order_shift_pixels(valid)-shift));
valid=valid&abs(diagnostic.order_shift_pixels-shift)<=max(.25,3*scatter);
diagnostic.accepted_orders=valid;
if nnz(valid)<3, return; end
shift=median(diagnostic.order_shift_pixels(valid));
scatter=1.4826*median(abs(diagnostic.order_shift_pixels(valid)-shift));
diagnostic.measured_shift_pixels=shift;
diagnostic.order_scatter_pixels=scatter;
if scatter>.5
    diagnostic.status='inconsistent_orders'; return
end
if abs(shift)>5
    diagnostic.status='offset_out_of_range'; return
end
for j=1:N
    field=['ofit' num2str(j)];
    if isfield(registered,field), registered.(field)=registered.(field)+shift; end
end
diagnostic.applied_shift_pixels=shift;
diagnostic.status='applied';
end
