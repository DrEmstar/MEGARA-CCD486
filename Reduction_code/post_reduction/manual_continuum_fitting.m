function [datacf]=manual_continuum_fitting(teff,logg,vsini,sigm,wstart,wstop,shift,firstorder,lastorder,data2,objectname)
%This program allows manual identification of the continuum level for each
%order. Try to keep the function smooth and fit the top sides of the line
%profiles.

%Initialise output
datacf=struct();

%Create synthetic spectrum
[synthwave,synthint,wavelengths,ele,elem,widths]=make_synth_spectrum(teff,logg,vsini,sigm,wstart,wstop);

%Apply velocity shift
synthwave=shift_velocity(synthwave,shift);

%Fit continuum order by order
for order_number=firstorder:lastorder
    intensity_field_name=['int' num2str(order_number)];
    wavelength_field_name=['wave' num2str(order_number)];
    continuum_intensity_field_name=['int' num2str(order_number)];
    continuum_wavelength_field_name=['w' num2str(order_number)];
    if ~isfield(data2,intensity_field_name) || ~isfield(data2,wavelength_field_name), warning(['Skipping order ' num2str(order_number) ' because data2.' intensity_field_name ' or data2.' wavelength_field_name ' is missing.']); continue; end
    intensity=data2.(intensity_field_name);
    wavelength=data2.(wavelength_field_name);
    temporary_synthetic_intensity=spline(synthwave,synthint,wavelength);
    if isvector(intensity)
        median_intensity=intensity;
    else
        median_intensity=median(intensity,1);
    end
    [selected_continuum_points,continuum_fit]=continuum_fitting_inline(wavelength,median_intensity,temporary_synthetic_intensity);
    datacf.(continuum_intensity_field_name)=continuum_fit;
    datacf.(continuum_wavelength_field_name)=wavelength;
end

%Save continuum fit
supplementary_directory=fullfile(pwd,'supplementary_data');
if ~isfolder(supplementary_directory), mkdir(supplementary_directory); end
custom_continuum_file_name=fullfile(supplementary_directory, ...
    ['datacf_' objectname '.mat']);
save(custom_continuum_file_name,'datacf')

end
