%MAKE_SYNSPEC_FILE Export a synthetic spectrum and its line list as text.
% Edit the atmosphere, broadening and wavelength settings below. This is a
% standalone export helper; it does not alter reduction products.

clearvars

sigm=0.05050;%FWHM full width at half maximum for Gaussian instrumental profile in Angstrom
wstart=3800;
wstop=8000;    
teff=5830;
logg=4.5;
vsini=-0.001;

[synthwave, synthint, wavelengths, ele, elem, widths] = make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop);

%ascii outputs

synthetic_spectrum=[synthwave synthint];
save('synthetic_spectrum_5830.dat','synthetic_spectrum','-ascii','-double')

outfile = 'line_list_5830.txt';
fid = fopen(outfile, 'w');
fprintf(fid, '# wavelength   width   ele   elem\n');
for i = 1:numel(wavelengths)
    fprintf(fid, '%.8f  %.6f  %s  %s\n', ...
        wavelengths(i), widths(i), ele{i}, elem{i});
end
fclose(fid);
