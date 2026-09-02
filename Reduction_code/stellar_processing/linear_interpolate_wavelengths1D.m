function wavelength=linear_interpolate_wavelengths1D(fwave1,fwave2,jd1,jd2,jdexp,m)

%Validate inputs
validate_wavelength_interpolation_inputs(fwave1,fwave2,jd1,jd2,jdexp,m);

%Initialise output
wavelength=struct();

%Calculate interpolation fraction
if jd2-jd1~=0
    P=(jdexp-jd1)/(jd2-jd1);
else
    P=0;
end

%Interpolate wavelength vectors order by order
for order_index=1:numel(m)
    order_number=m(order_index);
    order_field_name=['order_' num2str(order_number)];
    if ~isfield(fwave1,order_field_name), warning(['Skipping order ' num2str(order_number) ' because fwave1.' order_field_name ' is missing.']); continue; end
    if ~isfield(fwave2,order_field_name), warning(['Skipping order ' num2str(order_number) ' because fwave2.' order_field_name ' is missing.']); continue; end
    wave1=fwave1.(order_field_name);
    wave2=fwave2.(order_field_name);
    wave1=wave1(:);
    wave2=wave2(:);
    if numel(wave1)~=numel(wave2)
        warning(['Different number of wavelength elements in order ' num2str(order_number) ' - you should use the same flat-field for reduction to avoid this! Skipping this order.'])
        continue
    end
    wavelength.(order_field_name)=wave1+P*(wave2-wave1);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_wavelength_interpolation_inputs(fwave1,fwave2,jd1,jd2,jdexp,m)

%Check wavelength structures
if ~isstruct(fwave1), error('fwave1 must be a structure.'); end
if ~isstruct(fwave2), error('fwave2 must be a structure.'); end

%Check Julian dates
if ~isnumeric(jd1) || ~isscalar(jd1) || ~isfinite(jd1), error('jd1 must be a finite scalar numeric value.'); end
if ~isnumeric(jd2) || ~isscalar(jd2) || ~isfinite(jd2), error('jd2 must be a finite scalar numeric value.'); end
if ~isnumeric(jdexp) || ~isscalar(jdexp) || ~isfinite(jdexp), error('jdexp must be a finite scalar numeric value.'); end

%Check order list
if ~isnumeric(m) || ~isvector(m), error('m must be a numeric vector of order numbers.'); end
if any(~isfinite(m(:))) || any(m(:)~=round(m(:))), error('m must contain finite integer order numbers.'); end

end

