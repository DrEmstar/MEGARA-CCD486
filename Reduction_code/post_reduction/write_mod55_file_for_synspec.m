function write_mod55_file_for_synspec(wstart,wstop)
%change to appropriate synspec path

%Validate wavelength limits
if ~isnumeric(wstart) || ~isscalar(wstart) || ~isfinite(wstart), error('wstart must be a finite scalar numeric value.'); end
if ~isnumeric(wstop) || ~isscalar(wstop) || ~isfinite(wstop), error('wstop must be a finite scalar numeric value.'); end
if wstart>=wstop, error('wstart must be smaller than wstop.'); end

%Open SynSpec model-control file
file_identifier=fopen('mod.55','w+');
if file_identifier==-1, error('Could not open SynSpec input file: mod.55'); end
cleanup_object=onCleanup(@() fclose(file_identifier));

%Write SynSpec model-control settings
fprintf(file_identifier,'       0      55       0                    ! imode,idstd,iprint\n');
fprintf(file_identifier,'       0       0       0       0            ! imodel,ichang,interp,ichemc\n');
fprintf(file_identifier,'       0                                    ! Lyman quasi-mol\n');
fprintf(file_identifier,'       0       0       0       0       0\n');
fprintf(file_identifier,'     -23      24      26                    ! line broadening (H,HeI,HeII)\n');
fprintf(file_identifier,'    %.1f   %.1f    10       0  1.d-3    0.1\n',wstart,wstop);
fprintf(file_identifier,'      0\n');
fprintf(file_identifier,'      2.0\n');

end