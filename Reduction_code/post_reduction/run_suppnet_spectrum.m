function run_suppnet_spectrum(spectrum_files,varargin)
% RUN_SUPPNET_SPECTRUM  Call SUPPNET on one or more spectrum files via conda.
%
%   run_suppnet_spectrum("spec.txt")
%   run_suppnet_spectrum(["order76.txt","order77.txt"], 'weights','synth','sampling',[])
%
% specFiles : char, string, or cellstr; one or many paths to input spectra.

% -----------------------------
% Parse optional args
% -----------------------------

%Parse SUPPNET options
suppnet_options=parse_suppnet_options(varargin{:});

% -----------------------------
% Normalise specFiles -> string array
% -----------------------------

%Normalise input file list
spectrum_files=normalise_spectrum_file_list(spectrum_files);

if isempty(spectrum_files)
    error('No input spectrum files provided to run_suppnet_spectrum.');
end

%Check input files
for spectrum_file_index=1:numel(spectrum_files)
    if ~isfile(spectrum_files(spectrum_file_index)), error(['SUPPNET input spectrum file does not exist: ' char(spectrum_files(spectrum_file_index))]); end
end

% -----------------------------
% Paths to conda + SUPPNET script
% -----------------------------

%Check executable paths
if ~isfile(suppnet_options.suppnet_script), error(['SUPPNET script does not exist: ' char(suppnet_options.suppnet_script)]); end
if ~isfile(suppnet_options.conda_executable), error(['Conda executable does not exist: ' char(suppnet_options.conda_executable)]); end

% -----------------------------
% Build argument string
% -----------------------------

%Build sampling argument
if isempty(suppnet_options.sampling)
    sampling_argument='';
else
    sampling_argument=sprintf('--sampling %.6f',suppnet_options.sampling);
end

%Build quiet argument
if suppnet_options.quiet
    quiet_argument='--quiet';
else
    quiet_argument='';
end

%Build skip argument
if isempty(suppnet_options.skip) || suppnet_options.skip==0
    skip_argument='';
else
    skip_argument=sprintf('--skip %.6f',suppnet_options.skip);
end

%Build extra arguments
if strlength(string(suppnet_options.extra_arguments))==0
    extra_argument='';
else
    extra_argument=char(suppnet_options.extra_arguments);  % ensure char for sprintf
end

%Quote input files
quoted_spectrum_files=cell(1,numel(spectrum_files));
for spectrum_file_index=1:numel(spectrum_files)
    quoted_spectrum_files{spectrum_file_index}=shell_argument(spectrum_files(spectrum_file_index));
end
file_arguments=strjoin(quoted_spectrum_files,' ');

% Final command
command_text=sprintf('%s run -n %s python %s %s %s %s %s --weights %s %s',shell_argument(suppnet_options.conda_executable),shell_argument(suppnet_options.conda_environment),shell_argument(suppnet_options.suppnet_script),quiet_argument,sampling_argument,skip_argument,extra_argument,shell_argument(suppnet_options.weights),file_arguments);

%Run SUPPNET
fprintf('Running SUPPNET on %d file(s)...\n',numel(spectrum_files));
[command_status,command_output]=system(command_text);

%Check SUPPNET result
if command_status~=0
    error('SUPPNET failed with status %d:\n%s',command_status,command_output);
else
    fprintf('SUPPNET finished successfully.\n');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function suppnet_options=parse_suppnet_options(varargin)

%Accept structure or name-value options
if numel(varargin)==1 && isstruct(varargin{1})
    input_options=varargin{1};
    varargin=structure_to_name_value_arguments(input_options);
end

%Create parser
input_parser=inputParser;
addParameter(input_parser,'weights','synth');         % SUPPNET's internal weights: synth/active/emission
addParameter(input_parser,'sampling',0.1);            % wavelength sampling (or [] for native)
addParameter(input_parser,'conda_environment','suppnet-ktwo'); % conda env name
addParameter(input_parser,'env','');                  % legacy conda env name
addParameter(input_parser,'quiet',true);              % add --quiet flag
addParameter(input_parser,'skip',0);                  % SUPPNET --skip parameter (0 => omit)
addParameter(input_parser,'extra_arguments',"");      % extra CLI args as a string, e.g. "--something 3"
addParameter(input_parser,'extraArgs',"");            % legacy extra CLI args
addParameter(input_parser,'suppnet_script',"/full/path/to/suppnet/suppnet.py"); % configure during installation
addParameter(input_parser,'conda_executable',"/full/path/to/conda");             % configure during installation
parse(input_parser,varargin{:});

%Store parsed options
suppnet_options=input_parser.Results;

%Handle legacy option names
if strlength(string(suppnet_options.env))>0, suppnet_options.conda_environment=suppnet_options.env; end
if strlength(string(suppnet_options.extraArgs))>0, suppnet_options.extra_arguments=suppnet_options.extraArgs; end

%Validate option values
if isstring(suppnet_options.weights), suppnet_options.weights=char(suppnet_options.weights); end
if ~ischar(suppnet_options.weights) || ~any(strcmp(suppnet_options.weights,{'synth','active','emission'})), error('weights must be synth, active, or emission.'); end
if ~isempty(suppnet_options.sampling) && (~isnumeric(suppnet_options.sampling) || ~isscalar(suppnet_options.sampling) || ~isfinite(suppnet_options.sampling) || suppnet_options.sampling<=0), error('sampling must be empty or a finite positive scalar.'); end
if ~(islogical(suppnet_options.quiet) && isscalar(suppnet_options.quiet)) && ~(isnumeric(suppnet_options.quiet) && isscalar(suppnet_options.quiet) && any(suppnet_options.quiet==[0 1])), error('quiet must be true or false.'); end
suppnet_options.quiet=logical(suppnet_options.quiet);
if ~isempty(suppnet_options.skip) && (~isnumeric(suppnet_options.skip) || ~isscalar(suppnet_options.skip) || ~isfinite(suppnet_options.skip) || suppnet_options.skip<0), error('skip must be empty, zero, or a finite non-negative scalar.'); end

end

function name_value_arguments=structure_to_name_value_arguments(input_options)

%Convert structure to name-value arguments
field_names=fieldnames(input_options);
name_value_arguments=cell(1,2*numel(field_names));
for field_index=1:numel(field_names)
    name_value_arguments{2*field_index-1}=field_names{field_index};
    name_value_arguments{2*field_index}=input_options.(field_names{field_index});
end

end

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

end

function argument=shell_argument(value)

%Quote shell argument
if isnumeric(value), value=num2str(value); end
if isstring(value), value=char(value); end
if ~ischar(value), error('Shell argument must be a character vector, string scalar, or numeric scalar.'); end
value=strtrim(value);
if numel(value)>=2 && value(1)=='''' && value(end)=='''', argument=value; return; end
argument=['''' strrep(value,'''','''"''"''') ''''];

end
