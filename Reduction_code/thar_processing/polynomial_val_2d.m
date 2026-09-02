function z=polynomial_val_2d(c,x,y)
if nargin==0
    disp('usage: z = polynomial_val_2d(p,x,y)')
    z=[];
    return
end

% z = POLYVAL2D(V,sx,sy) Two dimensional polynomial evaluation
% If V is a matrix whose elements are the coefficients of a
% polynomial function of 2 variables, then POLYVAL2D(V,sx,sy)
% is the value of the polynomial evaluated at [sx,sy].  Row6X
% numbers in V correspond to powers of x, while column numbers
% in V correspond to powers of y. If sx and sy are matrices or
% vectors,the polynomial function is evaluated at all points
% in [sx, sy].
% If V is one dimensional, POLYVAL2D returns the same result as
% POLYVAL
% Use POLYFIT2D to generate appropriate polynomial matrices from
% f(x,y) data using a least squares method
% Perry W. Stout  June 28, 1995
% 4829 Rockland Way
% Fair Oaks, CA  95628
% (916) 966-0236
% Based on the Matlab function POLYVAL

% Polynomial evaluation c(x,y) is implemented using Horner's slick
% method.  Note use of the filter function to speed evaluation when
% the ordered pair [sx,sy] is single valued

%Set default y values
if nargin==2
    y=ones(size(x));
end

%Validate inputs
validate_polynomial_val_2d_inputs(c,x,y);

%Evaluate one-dimensional polynomial if needed
if isvector(c) && (size(c,1)==1 || size(c,2)==1)
    z=polyval(c,x);
    return
end

%Prepare coefficient dimensions
[number_x_coefficients,number_y_coefficients]=size(c);

%Evaluate polynomial using coefficient layout from polynomial_fit_2d
z=zeros(size(x));
for y_coefficient_index=1:number_y_coefficients
    x_coefficients=c(:,y_coefficient_index);
    x_polynomial_value=zeros(size(x));
    for x_coefficient_index=1:number_x_coefficients
        x_polynomial_value=x.*x_polynomial_value+x_coefficients(x_coefficient_index);
    end
    y_power=number_y_coefficients-y_coefficient_index;
    z=z+x_polynomial_value.*(y.^y_power);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_polynomial_val_2d_inputs(c,x,y)

%Check coefficient input
if ~isnumeric(c) || isempty(c), error('c must be a non-empty numeric coefficient array.'); end
if any(~isfinite(c(:))), error('c must contain only finite values.'); end

%Check coordinate inputs
if ~isnumeric(x) || ~isnumeric(y), error('x and y must be numeric.'); end
if ~isequal(size(x),size(y)), error('x and y must have the same dimensions.'); end
if any(~isfinite(x(:))) || any(~isfinite(y(:))), error('x and y must contain only finite values.'); end

end