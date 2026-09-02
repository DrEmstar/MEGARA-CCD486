function wave=make_wave_vector(thardata,wavfit,m)

%Validate inputs
validate_make_wave_vector_inputs(thardata,wavfit,m);

%Initialise output
wave=struct();

%Create wavelength vector for each order
for order_index=1:numel(m)
    order_number=m(order_index);
    order_field_name=['order_' num2str(order_number)];
    if ~isfield(thardata,order_field_name), warning(['Skipping missing ThAr order: thardata.' order_field_name]); continue; end
    if ~isfield(thardata.(order_field_name),'xax'), warning(['Skipping ThAr order ' num2str(order_number) ' because xax is missing.']); continue; end
    xpos=thardata.(order_field_name).xax;
    xpos=xpos(:);
    yy=(xpos-wavfit.ycen)/wavfit.yscale;
    xx=(order_number-wavfit.xcen)/wavfit.xscale;
    fitted_m_lambda=polynomial_val_2d(wavfit.p,repmat(xx,size(yy)),yy);
    fitted_m_lambda=fitted_m_lambda*wavfit.zscale+wavfit.zcen;
    wave.(order_field_name)=fitted_m_lambda/order_number;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_make_wave_vector_inputs(thardata,wavfit,m)

%Check ThAr data
if ~isstruct(thardata), error('thardata must be a structure.'); end

%Check wavelength fit
required_wavfit_fields={'p','xcen','ycen','zcen','xscale','yscale','zscale'};
for field_index=1:numel(required_wavfit_fields)
    if ~isfield(wavfit,required_wavfit_fields{field_index}), error(['wavfit is missing required field: ' required_wavfit_fields{field_index}]); end
end
if ~isnumeric(wavfit.p) || isempty(wavfit.p), error('wavfit.p must be a non-empty numeric coefficient matrix.'); end
if any(~isfinite(wavfit.p(:))), error('wavfit.p must contain only finite values.'); end
scale_values=[wavfit.xscale,wavfit.yscale,wavfit.zscale];
if any(~isfinite(scale_values)) || any(scale_values==0), error('wavfit scale values must be finite and non-zero.'); end
centre_values=[wavfit.xcen,wavfit.ycen,wavfit.zcen];
if any(~isfinite(centre_values)), error('wavfit centre values must be finite.'); end

%Check order numbers
if ~isnumeric(m) || ~isvector(m), error('m must be a numeric vector of order numbers.'); end
if any(~isfinite(m(:))) || any(m(:)==0), error('m must contain finite non-zero order numbers.'); end

end
