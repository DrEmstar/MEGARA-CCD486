

logg  = 4.0;
vsini = 50;
sigm  = 0.5;
wstart= 3600;
wstop = 8000;

cnt=0;
for teff = 7000 : 1000 : 20000
    cnt=cnt+1;
    eval(['[synthwave' num2str(cnt) ', synthint' num2str(cnt) ', wavelengths' num2str(cnt) ', ele' num2str(cnt) ', elem' num2str(cnt) ', widths' num2str(cnt) '] = make_synth(teff,logg,vsini,sigm,wstart,wstop);'])
end


teff = 7000 : 1000 : 20000;

rang = 4;
ratio= 0.1;
minim= 50;
for x=1:cnt
    fprintf('Good lines for Teff = %.0f\n\n',teff(x))
    eval(['find_good_lines_in_synthlist(wavelengths' num2str(x) ',widths' num2str(x) ',elem' num2str(x) ',rang,ratio,minim)'])
    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
end

synthetic_spectra = struct('synthwave1',synthwave1, 'synthint1',synthint1);


eval(['save(''synthspec_' objectname '.mat'',''synthetic_spectra'')'])


