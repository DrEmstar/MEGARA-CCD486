function info=display_directory_file_information_nooutput()
%This code queries the J*.fit files in a directory to dertermine their basic contents without a display

%Initialise output
info=struct();
info.list={};

%Find FITS files
file_list=dir('J*.fit');

%Open summary file
file_identifier=fopen('file_information.txt','w');
if file_identifier==-1, error('Could not open file_information.txt for writing.'); end
cleanup_file=onCleanup(@() fclose(file_identifier));

%Read file headers
for file_index=1:numel(file_list)
    file_root=file_list(file_index).name;
    file_root=file_root(1:end-4);
    if ~isvarname(file_root), error(['File root is not a valid MATLAB structure field name: ' file_root]); end
    header=get_header(file_list(file_index).name);
    temp=get_required_headers_from_header(header);
    temp=check_hercules_headers(temp);
    info.list{file_index}=file_root;
    info.(file_root)=temp;
    fprintf(file_identifier,'%s  %7s   %.4f    %d    %7.2f     %s\n',file_list(file_index).name,safe_string(temp.HERCEXPT),safe_number(temp.MJD_OBS),safe_integer(temp.HERCFIB),safe_number(temp.EXPTIME),safe_string(temp.OBJECT));
end

%Handle empty directories
if isempty(file_list)
    info.list=[];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output=safe_string(input_value)

%Format string value
if isempty(input_value)
    output='';
elseif isstring(input_value)
    output=char(input_value);
elseif ischar(input_value)
    output=input_value;
else
    output=char(string(input_value));
end

end

function output=safe_number(input_value)

%Format numeric value
if isempty(input_value) || ~isnumeric(input_value) || ~isscalar(input_value) || ~isfinite(input_value)
    output=NaN;
else
    output=input_value;
end

end

function output=safe_integer(input_value)

%Format integer value
if isempty(input_value) || ~isnumeric(input_value) || ~isscalar(input_value) || ~isfinite(input_value)
    output=NaN;
else
    output=round(input_value);
end

end

