function out=running_min(x,boxsize)

%Validate inputs
if ~isnumeric(x) || ~isvector(x), error('x must be a numeric vector.'); end
if ~isnumeric(boxsize) || ~isscalar(boxsize) || ~isfinite(boxsize) || boxsize<0, error('boxsize must be a finite scalar number >= 0.'); end

%Prepare window size
boxsize=round(boxsize);
window_size=2*boxsize+1;

%Calculate running minimum
out=movmin(x,window_size,'Endpoints','shrink');

end