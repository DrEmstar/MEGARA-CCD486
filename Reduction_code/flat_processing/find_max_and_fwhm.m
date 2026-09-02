function [centre,FWHM]=find_max_and_fwhm(xdata,ydata)

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
if numel(xdata)<3
    centre=0;
    FWHM=0;
    return
end

%Normalise profile
ydata=ydata-min(ydata);
maximum_ydata=max(ydata);
if ~isfinite(maximum_ydata) || maximum_ydata<=0
    centre=0;
    FWHM=0;
    return
end
ydata=ydata/maximum_ydata;

%Find expected centre
indcen=xdata(round(numel(xdata)/2));

%Find peaks
[pks,locs,widths,proms]=findpeaks(ydata,xdata,'MinPeakHeight',0.4,'WidthReference','halfheight');

%Select peak nearest expected centre
if numel(pks)>1
    nearest_peak_index=find_nearest(locs,indcen);
    centre=locs(nearest_peak_index);
    FWHM=widths(nearest_peak_index);
elseif numel(pks)==1
    centre=locs;
    FWHM=widths;
else
    centre=0;
    FWHM=0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nearest=find_nearest(vector,value)

%Find nearest index
[~,nearest]=min(abs(vector-value));

end