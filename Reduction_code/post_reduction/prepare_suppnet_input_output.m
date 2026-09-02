function [datacf]=prepare_suppnet_input_output(data2,objectname,firstorder,lastorder)
% PREPARE_SUPPNET_INPUT_OUTPUT
%   1. Checks if <objectname>_suppnet_datacf.mat already exists.
%      If so, loads it and returns immediately.
%   2. Otherwise:
%        - Builds median-per-order spectra
%        - Calls SUPPNET on them
%        - Loads .all outputs
%        - Constructs datacf structure (w76, int76, …)
%        - Saves <objectname>_suppnet_datacf.mat

%% ----------------------------------------------------
%  Check if datacf already exists
%% ----------------------------------------------------

%Set output paths
output_directory=fullfile(pwd,'suppnet_results');
output_datacf_file=fullfile(output_directory,sprintf('%s_suppnet_datacf.mat',objectname));

%Load cached continuum fit
if isfile(output_datacf_file)
    fprintf('Using existing suppnet fit\n');
    loaded_datacf_file=load(output_datacf_file);   % loads struct containing datacf
    if ~isfield(loaded_datacf_file,'datacf'), error(['Cached SUPPNET datacf file does not contain variable datacf: ' output_datacf_file]); end
    datacf=loaded_datacf_file.datacf;
    return
end

%% Ensure output directory exists

%Create output directory
if ~isfolder(output_directory)
    mkdir(output_directory);
end

%% ----------------------------------------------------
%  Build median spectra for SUPPNET
%% ----------------------------------------------------

%Initialise SUPPNET input file list
spectrum_files=strings(0,1);

%Write one median spectrum per order
for order_number=firstorder:lastorder
    wavelength_field_name=sprintf('wave%d',order_number);
    intensity_field_name=sprintf('int%d',order_number);
    if ~isfield(data2,wavelength_field_name) || ~isfield(data2,intensity_field_name), warning(['Skipping order ' num2str(order_number) ' because data2.' wavelength_field_name ' or data2.' intensity_field_name ' is missing.']); continue; end
    wavelength=data2.(wavelength_field_name);      % [N x 1]
    intensity=data2.(intensity_field_name);       % [N x M]
    if size(intensity,2)~=numel(wavelength), error(['data2.' intensity_field_name ' columns do not match data2.' wavelength_field_name ' length.']); end
    % Median across exposures (correct orientation)
    median_intensity=median(intensity,1,'omitnan')';   % [N x 1]
    % Output filename for SUPPNET input
    current_output_file=sprintf('%s_order_%d.txt',objectname,order_number);
    current_output_file=fullfile(output_directory,current_output_file);
    data_to_save=[wavelength(:),median_intensity(:)];
    writematrix(data_to_save,current_output_file,'Delimiter',' ','FileType','text');
    spectrum_files(end+1)=current_output_file;
end

%Stop if no spectra were written
if isempty(spectrum_files), error('No valid spectra were available for SUPPNET continuum fitting.'); end

%% ----------------------------------------------------
%  Run SUPPNET on all orders
%% ----------------------------------------------------

%Set SUPPNET options
suppnet_options=struct("quiet",true,"sampling",[],"weights","synth","skip",0,"extraArgs","");

%Run SUPPNET
run_suppnet_spectrum(spectrum_files,suppnet_options);

%% ----------------------------------------------------
%  Load SUPPNET .all outputs
%% ----------------------------------------------------

%Load SUPPNET results
normalised_data=load_suppnet_outputs(spectrum_files);

%Save SUPPNET results
output_results_file=fullfile(output_directory,sprintf('%s_suppnetresults.mat',objectname));
save(output_results_file,'normalised_data');

%% ----------------------------------------------------
%  Build datacf struct (w76, int76, etc.)
%% ----------------------------------------------------

%Initialise continuum-fit structure
datacf=struct();
orders=firstorder:lastorder;

%Convert SUPPNET output to datacf format
for order_index=1:numel(orders)
    order_number=orders(order_index);
    if order_index>numel(normalised_data), warning(['Missing SUPPNET output for order ' num2str(order_number) '.']); continue; end
    % Extract wavelength & normalized flux
    wavelength=normalised_data(order_index).wave;
    continuum_intensity=normalised_data(order_index).fit;   % SUPPNET normalized flux column
    % Field names
    wavelength_field_name=sprintf('w%d',order_number);
    continuum_intensity_field_name=sprintf('int%d',order_number);
    % Assign
    datacf.(wavelength_field_name)=wavelength(:)';
    datacf.(continuum_intensity_field_name)=continuum_intensity(:)';
end

%% ----------------------------------------------------
%  Save datacf
%% ----------------------------------------------------

%Save continuum-fit structure
output_datacf_file=fullfile(output_directory,sprintf('%s_suppnet_datacf.mat',objectname));
save(output_datacf_file,'datacf');

end