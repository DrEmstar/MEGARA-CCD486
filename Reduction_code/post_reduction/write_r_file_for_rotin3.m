function write_r_file_for_rotin3(vsini,sigm,wstart,wstop)
%change to appropriate synspec path

%Validate inputs
if ~isnumeric(vsini) || ~isscalar(vsini) || ~isfinite(vsini), error('vsini must be a finite scalar numeric value.'); end
if ~isnumeric(sigm) || ~isscalar(sigm) || ~isfinite(sigm), error('sigm must be a finite scalar numeric value.'); end
if ~isnumeric(wstart) || ~isscalar(wstart) || ~isfinite(wstart), error('wstart must be a finite scalar numeric value.'); end
if ~isnumeric(wstop) || ~isscalar(wstop) || ~isfinite(wstop), error('wstop must be a finite scalar numeric value.'); end
if vsini<0, error('vsini must be >= 0.'); end
if sigm<0, error('sigm must be >= 0.'); end
if wstart>=wstop, error('wstart must be smaller than wstop.'); end

%Open rotin3 input file
file_identifier=fopen('r.dat','w+');
if file_identifier==-1, error('Could not open rotin3 input file: r.dat'); end
cleanup_object=onCleanup(@() fclose(file_identifier));

%Write rotin3 settings
fprintf(file_identifier,' ''output.7''   ''output.17''    ''output.11'' \n');
%fprintf(fid,'     %.0f.    0.00       -1\n',vsini);high res
fprintf(file_identifier,'     %.0f.    0.05       0.1\n',vsini);
%fprintf(fid,'      %.2f             -1\n',sigm);high res
fprintf(file_identifier,'      %.2f             0\n',sigm);
fprintf(file_identifier,'    %.0f.   %.0f.      1\n',wstart,wstop);
fprintf(file_identifier,'\n');
fprintf(file_identifier,'\n');
fprintf(file_identifier,'\n');

end