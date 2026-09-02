function wavfit=wavelength_calibration(thfitinfo,order,air)

%Validate inputs
validate_wavelength_calibration_inputs(thfitinfo,order,air);

%Set constants
c=299792458; %speed of light (m/s)

%Build chosen line list
chosen.orders=order(thfitinfo(:,1));
chosen.mlam=order(thfitinfo(:,1)).*air(thfitinfo(:,1));
chosen.xpos=thfitinfo(:,4);
orig.orders=chosen.orders;
orig.mlam=chosen.mlam;
orig.xpos=chosen.xpos;

%Scale fit coordinates
xcen=min(chosen.orders);
ycen=min(chosen.xpos);
zcen=min(chosen.mlam);
xscale=max(abs(chosen.orders));
yscale=max(abs(chosen.xpos));
zscale=max(abs(chosen.mlam));
if xscale==0 || yscale==0 || zscale==0, error('Invalid wavelength calibration scaling value.'); end
x=(chosen.orders-xcen)/xscale;
y=(chosen.xpos-ycen)/yscale;
z=(chosen.mlam-zcen)/zscale;

%Fit 2D wavelength solution with iterative rejection
maximum_iterations=500;
minimum_lines_for_fit=50;
target_rms_metres_per_second=100;
for iteration_number=1:maximum_iterations
    if numel(z)<minimum_lines_for_fit, error('Too few ThAr lines remain for wavelength calibration.'); end
    p=polynomial_fit_2d(x,y,z,6,4);
    zfit=polynomial_val_2d(p,x,y);
    ZFIT=zfit*zscale+zcen;
    resid=chosen.mlam-ZFIT;
    er=std(resid,'omitnan');
    rver=er/mean(chosen.mlam,'omitnan')*c;
    if rver<=target_rms_metres_per_second
        fprintf('%.0f lines chosen, %.0f lines rejected\n\n',numel(ZFIT),numel(orig.mlam)-numel(ZFIT))
        wavfit.p=p;
        wavfit.xcen=xcen;
        wavfit.ycen=ycen;
        wavfit.zcen=zcen;
        wavfit.xscale=xscale;
        wavfit.yscale=yscale;
        wavfit.zscale=zscale;
        wavfit.rms_m_per_s=rver;
        wavfit.accepted_orders=chosen.orders;
        wavfit.accepted_xpos=chosen.xpos;
        wavfit.accepted_mlam=chosen.mlam;
        wavfit.original_orders=orig.orders;
        wavfit.original_xpos=orig.xpos;
        wavfit.original_mlam=orig.mlam;
        break
    end
    [~,worst_line_index]=max(abs(resid));
    chosen.xpos(worst_line_index)=[];
    chosen.orders(worst_line_index)=[];
    chosen.mlam(worst_line_index)=[];
    x(worst_line_index)=[];
    y(worst_line_index)=[];
    z(worst_line_index)=[];
end

%Stop if fit did not converge
if ~exist('wavfit','var'), error(['Wavelength calibration did not reach ' num2str(target_rms_metres_per_second) ' m/s after ' num2str(maximum_iterations) ' rejection iterations.']); end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_wavelength_calibration_inputs(thfitinfo,order,air)

%Check fitted line information
if ~isnumeric(thfitinfo) || size(thfitinfo,2)<4, error('thfitinfo must be a numeric matrix with at least four columns.'); end
if isempty(thfitinfo), error('thfitinfo is empty.'); end
if any(~isfinite(thfitinfo(:))), error('thfitinfo must contain only finite values.'); end

%Check reference line arrays
if ~isnumeric(order) || ~isvector(order), error('order must be a numeric vector.'); end
if ~isnumeric(air) || ~isvector(air), error('air must be a numeric vector.'); end
if numel(order)~=numel(air), error('order and air must have the same number of elements.'); end
if any(~isfinite(order(:))) || any(~isfinite(air(:))), error('order and air must contain only finite values.'); end

%Check line indices
line_indices=thfitinfo(:,1);
if any(line_indices<1) || any(line_indices>numel(order)) || any(line_indices~=round(line_indices)), error('thfitinfo(:,1) must contain valid integer indices into order and air.'); end

end