function [fwave1,fwave2]=make_fwave_variable(trueindbefore,trueindafter,fwave)
%this is required for parallel processing

%Validate inputs
if ~isstruct(fwave), error('fwave must be a structure.'); end
if isempty(trueindbefore) || isempty(trueindafter), error('trueindbefore and trueindafter must not be empty.'); end
if ~isnumeric(trueindbefore) || ~isnumeric(trueindafter), error('trueindbefore and trueindafter must be numeric indices.'); end

%Use first index if duplicate ThAr indices exist
before_index=trueindbefore(1);
after_index=trueindafter(1);

%Build field names
before_field_name=['th' num2str(before_index)];
after_field_name=['th' num2str(after_index)];

%Check wavelength solutions exist
if ~isfield(fwave,before_field_name), error(['fwave is missing field ' before_field_name '.']); end
if ~isfield(fwave,after_field_name), error(['fwave is missing field ' after_field_name '.']); end

%Return bracketing wavelength solutions
fwave1=fwave.(before_field_name);
fwave2=fwave.(after_field_name);

end