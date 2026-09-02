function [type,star]=identify_the_file(files,s)

%Validate inputs
if ~isstruct(files), error('files must be a structure.'); end
if ~isfield(files,'list'), error('files must contain a list field.'); end
if ~isnumeric(s) || ~isscalar(s) || s<1 || s>numel(files.list), error('s must be a valid scalar index into files.list.'); end

%Read file identifier
file_identifier=files.list{s};
if isstring(file_identifier), file_identifier=char(file_identifier); end
if ~isfield(files,file_identifier), error(['files does not contain metadata for ' file_identifier '.']); end

%Read exposure type
if isfield(files.(file_identifier),'HERCEXPT')
    type=files.(file_identifier).HERCEXPT; % either 'Thorium','White L' or 'Stellar'.
else
    type=[];
end

%Read object name
if isfield(files.(file_identifier),'OBJECT')
    star=files.(file_identifier).OBJECT; %'returns name of star for that file
else
    star=[];
end

end