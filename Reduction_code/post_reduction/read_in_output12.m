function [wavelengths,ele,elem,widths]=read_in_output12()
%Read SynSpec line-list output.

%Check output file exists
output_file_name='output.12';
if ~isfile(output_file_name), error('there must be a copy of the file output.12 in the current directory'); end

%Open output file
file_identifier=fopen(output_file_name,'r');
if file_identifier==-1, error(['Could not open SynSpec output file: ' output_file_name]); end
cleanup_object=onCleanup(@() fclose(file_identifier));

%Read line-list columns
line_list_columns=textscan(file_identifier,'%f %f %f %s %s %f %f %f %f %*s %f %f %f');

%Extract required columns
wavelengths=line_list_columns{3};
ele=line_list_columns{4};
ion=line_list_columns{5};
widths=line_list_columns{9};
elem=strcat(ele,ion);

%Remove zero-width lines
zero_width_lines=widths==0;
ele(zero_width_lines)=[];
elem(zero_width_lines)=[];
wavelengths(zero_width_lines)=[];
widths(zero_width_lines)=[];

end