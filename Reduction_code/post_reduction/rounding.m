function output=rounding(input,tonearest)

%Validate inputs
if ~isnumeric(input), error('input must be numeric.'); end
if ~isnumeric(tonearest) || ~isscalar(tonearest) || ~isfinite(tonearest) || tonearest==0, error('tonearest must be a finite non-zero scalar numeric value.'); end

%Round input to nearest requested value
output=round(input/tonearest)*tonearest;

end
