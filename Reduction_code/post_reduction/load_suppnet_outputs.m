function normalised_data=load_suppnet_outputs(spectrum_files)
% LOAD_SUPPNET_OUTPUTS  Load normalised spectra from SUPPNET .all files.
%
%   normData = load_suppnet_outputs(specFiles)
%
% specFiles : char, string, or cellstr
%             The *input* spectrum files you passed to run_suppnet_spectrum.
%
% Returns:
%   normData is a struct array with fields:
%       .inputFile   - original input file path (the .txt)
%       .outputFile  - .all file path
%       .wave        - wavelength array       (column 1)
%       .flux        - normalised flux array  (column 2, adjust if needed)

% Normalise input to string array
spectrum_files=normalise_spectrum_file_list(spectrum_files);

%Initialise output structure
number_of_spectrum_files=numel(spectrum_files);
normalised_data=repmat(struct('inputFile',[],'outputFile',[],'wave',[],'influx',[],'outflux',[],'fit',[],'flux',[]),number_of_spectrum_files,1);

%Load each SUPPNET output file
for spectrum_file_index=1:number_of_spectrum_files
    input_file=spectrum_files(spectrum_file_index);
    [output_directory,base_file_name,~]=fileparts(input_file);

    % SUPPNET output file: same base, .all extension
    output_file=fullfile(output_directory,base_file_name+".all");

    %Check output file exists
    if ~isfile(output_file)
        error('Could not find SUPPNET output for "%s". Expected "%s".',input_file,output_file);
    end

    % Tell MATLAB it's a text file regardless of extension
    output_data=readmatrix(output_file,'FileType','text');

    %Check output data
    if isempty(output_data)
        error('Output file "%s" is empty.',output_file);
    end
    if size(output_data,2)<5
        error('Output file "%s" must have at least five columns because this code expects wave=1, influx=2, outflux=3, fit=5.',output_file);
    end

    %Store output data
    normalised_data(spectrum_file_index).inputFile=char(input_file);
    normalised_data(spectrum_file_index).outputFile=char(output_file);
    normalised_data(spectrum_file_index).wave=output_data(:,1);
    normalised_data(spectrum_file_index).influx=output_data(:,2);
    normalised_data(spectrum_file_index).outflux=output_data(:,3);
    normalised_data(spectrum_file_index).fit=output_data(:,5);
    normalised_data(spectrum_file_index).flux=output_data(:,5);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spectrum_files=normalise_spectrum_file_list(spectrum_files)

%Normalise spectrum file input
if ischar(spectrum_files) || isstring(spectrum_files)
    spectrum_files=string(spectrum_files);
elseif iscellstr(spectrum_files)
    spectrum_files=string(spectrum_files);
else
    error('specFiles must be char, string, or cellstr.');
end
spectrum_files=spectrum_files(:);

%Check input files
if isempty(spectrum_files)
    error('No input spectrum files provided to load_suppnet_outputs.');
end

end