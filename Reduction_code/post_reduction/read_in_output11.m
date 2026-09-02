function [synthwave,synthint]=read_in_output11()
%Read SynSpec synthetic spectrum output.

%Check output file exists
output_file_name='output.11';
if ~isfile(output_file_name), error('there must be a copy of the file output.11 in the current directory'); end

%Open output file
file_identifier=fopen(output_file_name,'r');
if file_identifier==-1, error(['Could not open SynSpec output file: ' output_file_name]); end
cleanup_object=onCleanup(@() fclose(file_identifier));

%Read wavelength and intensity columns
spectrum_columns=textscan(file_identifier,'%f %f');

%Extract spectrum
synthwave=spectrum_columns{1};
synthint=spectrum_columns{2};

%Check output was read
if isempty(synthwave) || isempty(synthint), error(['No spectrum data were read from ' output_file_name]); end

end