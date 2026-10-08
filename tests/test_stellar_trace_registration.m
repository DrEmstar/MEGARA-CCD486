function test_stellar_trace_registration
% Deterministic centroid, outlier, low-signal and detector-boundary tests.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','Reduction_code')));
[a,s,b,f]=profiles(2.4);
[r,d]=register_stellar_trace(a,s,b,f);
assert(strcmp(d.status,'applied'));assert(abs(d.applied_shift_pixels-2.4)<.01);
assert(max(abs(r.ofit1-a.ofit1-2.4))<.01);
% Master profile is asymmetric relative to nominal trace; compare to it.
[a,s,b,f]=profiles(-1.2);[~,d]=register_stellar_trace(a,s,b,f);
assert(abs(d.applied_shift_pixels+1.2)<.01);
% One contaminated order must not drive the global shift.
s.order_5.data=fliplr(s.order_5.data);[~,d]=register_stellar_trace(a,s,b,f);
assert(abs(d.applied_shift_pixels+1.2)<.01);assert(~d.accepted_orders(5));
% Disabled mode preserves the trace exactly.
a.stellar_trace_registration_enabled=false;[r,d]=register_stellar_trace(a,s,b,f);
assert(isequal(r,a)&&strcmp(d.status,'disabled'));
a=rmfield(a,'stellar_trace_registration_enabled');
for j=1:5,s.(['order_' num2str(j)]).data(:)=0;end
[r,d]=register_stellar_trace(a,s,b,f);assert(isequal(r,a)&&d.applied_shift_pixels==0);
[a,s,b,f]=profiles(6);[r,d]=register_stellar_trace(a,s,b,f);
assert(isequal(r,a)&&strcmp(d.status,'offset_out_of_range'));
[a,s,b,f]=profiles(1);
for j=1:4,f=rmfield(f,['order_' num2str(j)]);end
[r,d]=register_stellar_trace(a,s,b,f);assert(isequal(r,a)&&strcmp(d.status,'insufficient_orders'));
% Large disagreement between orders rejects a global translation.
[a,s,b,f]=profiles(0);
for j=1:5
 offs=-15:15;s.(['order_' num2str(j)]).data=repmat(1000*exp(-.5*((offs-.8-(j-3))/2).^2),200,1);
end
[r,d]=register_stellar_trace(a,s,b,f);assert(isequal(r,a)&&strcmp(d.status,'inconsistent_orders'));
% Nonfinite columns and mismatched calibration geometry are rejected.
[a,s,b,f]=profiles(.7);s.order_1.data(:)=NaN;f.order_2.xax=f.order_2.xax+1;
[~,d]=register_stellar_trace(a,s,b,f);assert(nnz(d.accepted_orders)==3);assert(abs(d.applied_shift_pixels-.7)<.01);
% End-to-end synthetic extraction, including cosmic filtering and weights geometry.
[a,~,~,f]=profiles(1);a.WIDTHFACTOR=3;
image=zeros(230,200);yp=(1:230)';
for j=1:5
 a.(['width' num2str(j)])=5;a.(['points' num2str(j)])=(1:200)';
 c=30+40*(j-1);a.(['ofit' num2str(j)])=repmat(c,200,1);
 image=image+1000*exp(-.5*((yp-c-1.8)/2).^2)*ones(1,200);
end
[star,back,d]=extract_registered_stellar_orders(a,image,zeros(size(image)),f,0);
assert(strcmp(d.status,'applied'));assert(abs(d.applied_shift_pixels-1)<.02);
assert(isempty(d.fallback_orders));assert(isequal(star.order_1.xax,f.order_1.xax));
assert(all(back.order_1.data(:)==0));
expected=0;for j=1:5,expected=expected+nnz(star.(['order_' num2str(j)]).cosmics);end
assert(star.cosmic_count==expected);
a.stellar_trace_registration_enabled=false;
[fixed,~,d]=extract_registered_stellar_orders(a,image,zeros(size(image)),f,0);
original=extract_all_orders_no_background(a,image);original.numords=a.numords;original=remove_cosmics(original,0);
assert(isequaln(fixed,original)&&strcmp(d.status,'disabled'));
a=rmfield(a,'stellar_trace_registration_enabled');
% Edge clipping changes xax: keep original aperture rather than misindex flat/wavelength.
a.ofit1(:)=16;image=zeros(230,200);
for j=1:5
 c=a.(['ofit' num2str(j)])(1);image=image+1000*exp(-.5*((yp-c+.2)/2).^2)*ones(1,200);
end
[star,~,d]=extract_registered_stellar_orders(a,image,zeros(size(image)),f,0);
assert(d.applied_shift_pixels<-.9);assert(ismember(1,d.fallback_orders));assert(numel(star.order_1.xax)==200);
disp('All stellar trace registration tests passed.');
end
function [a,s,b,f]=profiles(shift)
a=struct('numords',5);s=struct;b=struct;f=struct;offs=-15:15;
for j=1:5
 c=30+40*(j-1);a.(['ofit' num2str(j)])=repmat(c,200,1);
 q=struct('xax',(1:200)','ypositions',repmat(c,200,1),'yax',repmat(c+offs,200,1));
 q.data=repmat(1000*exp(-.5*((offs-.8)/2).^2),200,1);q.extraction_inds=true(1,31);f.(['order_' num2str(j)])=q;
 q.data=repmat(1000*exp(-.5*((offs-.8-shift)/2).^2),200,1);s.(['order_' num2str(j)])=q;
 q.data=zeros(200,31);b.(['order_' num2str(j)])=q;
end
end
