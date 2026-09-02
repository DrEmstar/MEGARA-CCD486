function [x,resnorm,residual,exitflag,fitted]=fit_th_line(xdata,ydata,lb,ub,startpoint)

%Validate inputs
validate_thorium_line_fit_inputs(xdata,ydata,lb,ub,startpoint);

%Prepare vectors
xdata=xdata(:);
ydata=ydata(:);

%fit simple Gaussian to data
fun=@(x,xdata) x(1)*exp(-(((xdata-x(2))*(2*sqrt(log(2)))/max(abs(x(3)),eps)).^2))+x(4)+x(5)*xdata;%,'x','xdata';

%Fit Gaussian model
options=optimset('Display','off');
try
    [x,resnorm,residual,exitflag]=lsqcurvefit(fun,startpoint,xdata,ydata,lb,ub,options);
    fitted=fun(x,xdata);
catch fit_error
    warning(['ThAr line fit failed: ' fit_error.message])
    x=nan(size(startpoint));
    resnorm=Inf;
    residual=nan(size(ydata));
    exitflag=-999;
    fitted=nan(size(ydata));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_thorium_line_fit_inputs(xdata,ydata,lb,ub,startpoint)

%Check data vectors
if ~isnumeric(xdata) || ~isvector(xdata), error('xdata must be a numeric vector.'); end
if ~isnumeric(ydata) || ~isvector(ydata), error('ydata must be a numeric vector.'); end
if numel(xdata)~=numel(ydata), error('xdata and ydata must have the same number of elements.'); end
if any(~isfinite(xdata(:))) || any(~isfinite(ydata(:))), error('xdata and ydata must contain only finite values.'); end
if numel(xdata)<5, error('At least five data points are required to fit a ThAr line.'); end

%Check bounds and starting point
if ~isnumeric(lb) || ~isnumeric(ub) || ~isnumeric(startpoint), error('lb, ub, and startpoint must be numeric vectors.'); end
if numel(lb)~=5 || numel(ub)~=5 || numel(startpoint)~=5, error('lb, ub, and startpoint must each contain five parameters.'); end
if any(~isfinite(lb(:))) || any(~isfinite(ub(:))) || any(~isfinite(startpoint(:))), error('lb, ub, and startpoint must contain only finite values.'); end
if any(lb(:)>ub(:)), error('Each lower bound must be <= the corresponding upper bound.'); end
if any(startpoint(:)<lb(:)) || any(startpoint(:)>ub(:)), error('startpoint must lie within lb and ub.'); end

end