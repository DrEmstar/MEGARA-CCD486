function [temp,jdexp]=make_temp_variable(info,s)
%this is required by the parallel loop

%Validate inputs
if ~isstruct(info), error('info must be a structure.'); end
if ~isfield(info,'list'), error('info must contain a list field.'); end
if ~isnumeric(s) || ~isscalar(s) || ~isfinite(s) || s<1 || s>numel(info.list), error('s must be a valid scalar index into info.list.'); end

%Read file identifier
file_identifier=info.list{s};
if isstring(file_identifier), file_identifier=char(file_identifier); end
if ~isfield(info,file_identifier), error(['info does not contain metadata for ' file_identifier '.']); end

%Read exposure type
if isfield(info.(file_identifier),'HERCEXPT')
    temp=info.(file_identifier).HERCEXPT;
else
    temp=[];
end

%Read Julian date
if isfield(info.(file_identifier),'MJD_OBS')
    jdexp=info.(file_identifier).MJD_OBS;
else
    jdexp=[];
end

end