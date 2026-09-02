function ldelfn_array=logdelfn_2_logdelfn_array(ldelfn,logwavelength)
%Convert logarithmic delta-function line list to a sampled array.

%Validate inputs
validate_log_delta_function_inputs(ldelfn,logwavelength);

%Initialise sampled line-list array
ldelfn_array=zeros(size(logwavelength));

%Build nearest-neighbour bin edges
logwavelength_vector=logwavelength(:);
bin_edges=[-Inf;(logwavelength_vector(1:end-1)+logwavelength_vector(2:end))/2;Inf];

%Assign lines to nearest wavelength pixel
nearest_indices=discretize(ldelfn(:,1),bin_edges);

%Sum line weights in each pixel
summed_weights=accumarray(nearest_indices,ldelfn(:,2),[numel(logwavelength_vector) 1],@sum,0);

%Restore original wavelength-axis shape
ldelfn_array(:)=summed_weights;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_log_delta_function_inputs(ldelfn,logwavelength)

%Check line-list input
if ~isnumeric(ldelfn) || size(ldelfn,2)~=2, error('ldelfn must be a numeric array with two columns: log wavelength and weight.'); end
if isempty(ldelfn), error('ldelfn is empty.'); end
if any(~isfinite(ldelfn(:))), error('ldelfn must contain only finite values.'); end

%Check wavelength grid
if ~isnumeric(logwavelength) || ~isvector(logwavelength), error('logwavelength must be a numeric vector.'); end
if isempty(logwavelength), error('logwavelength is empty.'); end
if any(~isfinite(logwavelength(:))), error('logwavelength must contain only finite values.'); end
if any(diff(logwavelength(:))<=0), error('logwavelength must be strictly increasing.'); end

end