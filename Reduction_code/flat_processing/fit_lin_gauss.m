function [centre,widths,sse,F]=fit_lin_gauss(xdata,ydata,start_point)

%this function is based on the "Fitting a Curve to Data" help topic in the
%MATLAB help

%Validate inputs
if ~isnumeric(xdata) || ~isvector(xdata), error('xdata must be a numeric vector.'); end
if ~isnumeric(ydata) || ~isvector(ydata), error('ydata must be a numeric vector.'); end
if numel(xdata)~=numel(ydata), error('xdata and ydata must have the same number of elements.'); end

%Prepare vectors
xdata=xdata(:);
ydata=ydata(:);
valid_values=isfinite(xdata) & isfinite(ydata);
xdata=xdata(valid_values);
ydata=ydata(valid_values);

%Return failed fit if too few points remain
if numel(xdata)<5
    centre=NaN;
    widths=NaN;
    sse=Inf;
    F=NaN(size(xdata));
    return
end

% Call fminsearch with a fixed starting point.
% this reduction version always has a specified width startpoint ...
if ~exist('start_point','var') || isempty(start_point)
    [maximum_ydata,maximum_index]=max(ydata);
    start_point=[maximum_ydata,xdata(maximum_index),0,min(ydata),5]; %assuming centred and scaled then these should be reasonable starting points
end

%parameters are:
%[amp xshift slope yshift width];

%Fit linear Gaussian model
model=@lingaussfun;
options=optimset('Display','off');
estimates=fminsearch(model,start_point,options);

%Evaluate fitted model
F=evaluate_linear_gaussian(estimates,xdata);

%Store fitted parameters
centre=estimates(2);
widths=abs(estimates(5));
sse=lingaussfun(estimates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitted_curve=evaluate_linear_gaussian(parameters,x_values)

%Evaluate model curve
amplitude=parameters(1);
xshift=parameters(2);
slope=parameters(3);
yshift=parameters(4);
width=abs(parameters(5));
if width==0 || ~isfinite(width), fitted_curve=Inf(size(x_values)); return; end
fitted_curve=amplitude*exp(-(((x_values-xshift)*2*sqrt(log(2))/width).^2))+slope*x_values+yshift;

end

function sse=lingaussfun(params)%[sse, FittedCurve] = lingaussfun(params)

%Calculate sum of squared errors
if numel(params)<5 || any(~isfinite(params)) || params(5)==0
    sse=Inf;
    return
end
a=params(1);
b=params(2);
c=params(3);
d=params(4);
e=abs(params(5));
FittedCurve=a*exp(-(((xdata-b)*2*sqrt(log(2))/e).^2))+c*xdata+d; %the best fit curve shape - in this case a gaussian on a sloped line
ErrorVector=FittedCurve-ydata; %least-squares error
sse=sum(ErrorVector.^2); %used for the fminsearch optimisation

end

end