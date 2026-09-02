function write_new_file_for_synspec(teff,logg)
%change to appropriate synspec path

%Store current directory
starting_directory=pwd;
cleanup_object_directory=onCleanup(@() cd(starting_directory));

%Move to SynSpec directory
synspec_pointer_file=which('synspec_pointer.m');
if isempty(synspec_pointer_file), error('Could not find synspec_pointer.m on the MATLAB path.'); end
synspec_directory=fileparts(synspec_pointer_file);
cd(synspec_directory)

%Validate inputs
if ~isnumeric(teff) || ~isscalar(teff) || ~isfinite(teff), error('teff must be a finite scalar numeric value.'); end
if ~isnumeric(logg) || ~isscalar(logg) || ~isfinite(logg), error('logg must be a finite scalar numeric value.'); end

%Write SynSpec input file
file_identifier=fopen('new','w+');
if file_identifier==-1, error('Could not open SynSpec input file: new'); end
cleanup_object_file=onCleanup(@() fclose(file_identifier));

%Write stellar parameters
if teff<10000
    fprintf(file_identifier,'       %.0f.          %.2f',teff,logg);
else
    fprintf(file_identifier,'      %.0f.          %.2f',teff,logg);
end

end