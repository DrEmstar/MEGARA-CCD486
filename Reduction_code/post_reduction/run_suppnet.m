function run_suppnet(inFile, outFile, varargin)
% RUN_SUPPNET  Call SUPPNET from MATLAB via your conda env.
%
% Example:
%   run_suppnet('data/input.fasta','results/out.txt','--model','configs/model.yaml');

    % Make sure required paths are available to the subprocess
    addpaths = ':/Users/ejb555/bin:/opt/anaconda3/bin';
    setenv('PATH', [getenv('PATH') addpaths]);

    % Quote paths in case they contain spaces
    inFile  = string(inFile);
    outFile = string(outFile);
    q = @(s) sprintf('"%s"', s);

    % Build CLI (edit subcommand/flags per SUPPNET usage)
    args = strjoin(string(varargin), ' ');
    cli  = sprintf('/Users/ejb555/bin/SUPPNET --path %s %s', q(inFile), args);


    % Prefer conda run to ensure the right env
    cmd = sprintf('conda run -n suppnet-env %s', cli);

    % Execute
    [status, cmdout] = system(cmd);

    % Surface output and errors
    fprintf('%s\n', cmdout);
    if status ~= 0
        error('SUPPNET failed (status %d). See output above.', status);
    end
end
