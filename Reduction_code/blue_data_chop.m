function data=blue_data_chop(data,value)
%chop off blue orders (first 'value' pixel rows)

%Validate inputs
if ~isnumeric(data) || ndims(data)~=2, error('data must be a numeric 2D array.'); end
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value<0, error('value must be a finite scalar number >= 0.'); end

%Round chop value
value=round(value);

%Return unchanged data if no chop requested
if value==0
    return
end

%Check chop size
if value>=size(data,1), error('blue_data_chop value removes all rows from data. Reduce blue_data_chop_value.'); end

%chop off blue orders (first 'value' pixel rows)
data(1:value,:)=[];

end
