function p=polynomial_fit_2d(x,y,z,n,m)
if nargin==0
    disp('usage: p = polynomial_fit_2d(x,y,z,n,m)')
    p=[];
    return
end

% P= POLYFIT2D(x,y,z,n,m) finds the coefficients of a 
%  polynomial function of 2 variables formed from the
%  data in vectors x and y of degrees n and m, respectively,
%  that fit	the data in vector z in a least-squares sense.  
%
% The regression problem is formulated in matrix format as:
%
%	z = A*P   So that if the polynomial is cubic in  
%             x and linear in y, the problem becomes:
%
%	z = [y.*x.^3  y.*x.^2  y.*x  y x.^3  x.^2  x  ones(length(x),1)]*
%	    [p31 p21 p11 p01 p30 p20 p10 p00]'                      
%		
%  Note that the various xy products are column vectors of length(x).
%
%  The coefficents of the output p    
%  matrix are arranged as shown:
%
%      p31 p30 
%      p21 p20 
%      p11 p10 
%      p01 p00
%
% The indices on the elements of p correspond to the 
% order of x and y associated with that element.
%
% For a solution to exist, the number of ordered 
% triples [x,y,z] must equal or exceed (n+1)*(m+1).
% Note that m or n may be zero.
%
% To evaluate the resulting polynominal function,
% use POLYVAL2D.
%
% Perry W. Stout  June 29, 1995
% 4829 Rockland Way
% Fair Oaks, CA  95628
% (916) 966-0236
% Based on the Matlab function polyfit.

%Validate inputs
validate_polynomial_fit_2d_inputs(x,y,z,n,m);

%Convert inputs to columns
x=x(:);
y=y(:);
z=z(:);  % Switches vectors to columns--matrices, too

%Check number of data points
number_coefficients=(n+1)*(m+1);
if numel(x)<number_coefficients
    error('Number of points must equal or exceed order of polynomial function.')
end

%Set polynomial dimensions
x_degree_count=n+1;
y_degree_count=m+1; % Increments n and m to equal row, col numbers of p.

% Construct the extended Vandermonde matrix, containing all xy products.
design_matrix=zeros(numel(x),x_degree_count*y_degree_count);
for y_power_index=1:y_degree_count
    for x_power_index=1:x_degree_count
        design_matrix(:,x_power_index+(y_power_index-1)*x_degree_count)=(x.^(x_degree_count-x_power_index)).*(y.^(y_degree_count-y_power_index));
    end
end

%Solve least-squares system
coefficient_vector=design_matrix\z;

% Reform p as a matrix.
p=reshape(coefficient_vector,x_degree_count,y_degree_count);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_polynomial_fit_2d_inputs(x,y,z,n,m)

%Check coordinate vectors
if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(z), error('x, y, and z must be numeric.'); end
if numel(x)~=numel(y) || numel(z)~=numel(y), error('X, Y,and Z vectors must be the same size'); end
if any(~isfinite(x(:))) || any(~isfinite(y(:))) || any(~isfinite(z(:))), error('x, y, and z must contain only finite values.'); end

%Check polynomial degrees
if ~isnumeric(n) || ~isscalar(n) || ~isfinite(n) || n<0 || n~=round(n), error('n must be a non-negative integer scalar.'); end
if ~isnumeric(m) || ~isscalar(m) || ~isfinite(m) || m<0 || m~=round(m), error('m must be a non-negative integer scalar.'); end

end