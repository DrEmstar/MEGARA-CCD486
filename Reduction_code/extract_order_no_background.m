function extracted=extract_order_no_background(points,ofit,width,img,widthfactor)

% Fractional-trace extraction using shape-preserving interpolation.
% Uses pchip instead of linear interpolation to reduce sub-pixel phase artefacts.

%Validate inputs
validate_extract_order_inputs(points,ofit,width,img,widthfactor);

%Prepare image geometry
szimg=size(img);
xax=round(points(:));

% keep fitted trace centre fractional
ypositions=ofit(:);

%Build extraction offsets
offsets=ceil(-width*widthfactor):floor(width*widthfactor);
nrows=numel(offsets);
if nrows==0, error('Extraction aperture is empty. Check width and widthfactor.'); end

%Preallocate extraction arrays
sxax=zeros(size(xax));
syax=zeros(numel(xax),nrows);
sdata=zeros(numel(xax),nrows);
sypositions=zeros(size(xax));

%Build detector row axis
ypix=(1:szimg(1))';

%Extract valid columns
cnt=0;
for xpos=1:numel(xax)
    xc=xax(xpos);
    yc=ypositions(xpos);
    if xc<1 || xc>szimg(2) || ~isfinite(yc), continue; end
    ygrid=yc+offsets;
    if min(ygrid)>=1 && max(ygrid)<=szimg(1)
        cnt=cnt+1;
        sxax(cnt)=xc;
        syax(cnt,:)=ygrid;
        % Changed from 'linear' to 'pchip'
        sdata(cnt,:)=interp1(ypix,double(img(:,xc)),ygrid,'pchip');
        sypositions(cnt)=yc;
    end
end

%Warn if extraction failed
if cnt==0
    warning('No valid pixels were extracted for this order.')
end

%Store extracted order
extracted.xax=sxax(1:cnt);
extracted.yax=syax(1:cnt,:);
extracted.data=sdata(1:cnt,:);
extracted.ypositions=sypositions(1:cnt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_extract_order_inputs(points,ofit,width,img,widthfactor)

%Check trace inputs
if ~isnumeric(points) || ~isvector(points), error('points must be a numeric vector.'); end
if ~isnumeric(ofit) || ~isvector(ofit), error('ofit must be a numeric vector.'); end
if numel(points)~=numel(ofit), error('points and ofit must have the same number of elements.'); end
if any(~isfinite(points(:))) || any(~isfinite(ofit(:))), error('points and ofit must contain only finite values.'); end

%Check extraction geometry
if ~isnumeric(width) || ~isscalar(width) || ~isfinite(width) || width<=0, error('width must be a finite positive scalar.'); end
if ~isnumeric(widthfactor) || ~isscalar(widthfactor) || ~isfinite(widthfactor) || widthfactor<=0, error('widthfactor must be a finite positive scalar.'); end

%Check image
if ~isnumeric(img) || ndims(img)~=2, error('img must be a numeric 2D array.'); end
if any(~isfinite(img(:))), warning('img contains non-finite values.'); end

end