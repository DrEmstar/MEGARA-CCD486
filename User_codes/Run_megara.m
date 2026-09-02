%Welcome to MEGARA. This is Version 2.1

clearvars
clc

%USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

options=struct();

options.objectname=['hd_48501'];

%If objectname is empty then all stellar images are processed
%to run all MARSDEN EXOCOMET stars use 'exocomet'
%to run all MUSICIAN stars use 'musician'
%to run all MUSICIAN TESS stars use 'tess'
%to run a list of targets define a structure array, e.g. objectname.target1='hd_000001';objectname.target2='hd_000002'; ...

%Months/Runs to reduce
options.reduction_folders={['tutorial_data']};
%A colon-separated cell array of the directory names to be reduced, e.g. reduction_folders={'Jan 2018'; 'Feb 2018'};
%If reduction_folders is empty then all folders in the directory where pointer_raw_data_directory.m is located are reduced

%Input stellar values (these can be left blank for MUSICIAN or EXOCOMET targets)
options.teff=[]; %Effective temperature in Kelvin
options.logg=[]; %Surface gravity as log10(g) in cgs units
options.vsini=[]; %Projected rotational velocity in km/s

% options.teff=8000; %Effective temperature in Kelvin
% options.logg=4.0; %Surface gravity as log10(g) in cgs units
% options.vsini=60; %Projected rotational velocity in km/s

%USER OPTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose true for on and false for off.

options.star_census=true; %Do you want to redo census of observations? (default true)
options.skip_reduction=false %Do you want to skip reduction? (default false)
options.skip_post_reduction=false; %Do you want to skip processing and post-reduction? (default false)
options.skip_manual_file_check=false; %Skip the Good/Bad interface and retain ambiguous files for unattended reductions (default false)

options.chunked_post_reduction=false; %Process already-reduced spectra internally in chunks. Recommended for large data sets.
options.post_reduction_chunk_size=500; %Number of already-reduced spectra per post-reduction chunk.
options.extended_wavelength_range=false; %Default false gives 3869-6866 Angstrom. Set true for approximately 3750-8950 Angstrom; manual review and additional processing may be required at the extreme blue and red ends.

options.blue_data_chop_value=300; %Number of pixels to cut off the blue end of the spectrum due to incomplete flats (default=500, 300 for Ca H&K)
options.apply_barycentric_correction=true; %Do you want to apply the barycentric correction? (default true)

options.apply_continuum_fit=true; %Apply a continuum fit to normalise the data (default=true)
%Continuum fitting uses a previous manual fit in the reduced data directory for the star by default. If this does not exist then one of the options below must be set.
options.manual_continuum_fit=false; %Do you want to make a new manual continuum fit? Use this for precision work. This overwrites previous manual fits; save those elsewhere first.
options.template_continuum_fit=false; %Apply the interpolated SpT continuum grid directly before order merging. Requires teff and the standard wavelength range; SUPPNet is not applied.
options.automatic_continuum_fit=true; %Use the two-stage automatic fit: synthetic-guided order normalisation, harmonised merging, then SUPPNet on the full merged spectrum.

options.order_merge=true; %Do you want to merge orders? (default true) Requires apply_continuum_fit=true

options.overwrite_full_data=true; %Use the existing order-merged spectra while regenerating the LSD profiles (default true)
options.overwrite_post_reduction=true; %Do you want to redo regularised LSD profiles and related measurements? (default true)

options.calculate_lsd_profile=true; %Do you want to calculate regularised LSD profiles? Requires order_merge=true. This overwrites any previous regularised LSD profile for this object.
options.radial_velocity_measurement=true; %Do you want to measure radial velocity and vsini from the regularised LSD profiles? Requires calculate_lsd_profile=true.
options.extend_velocities=false; %Extends the LSD velocity axes from +/- 200 km/s to +/- 400 km/s (default false)

options.suppress_figures=false; %Do you want to suppress figures during running of MEGARA? (default false)

options.ascii_output=false; %Do you want your output files in ASCII format?
options.fits_output=false; %Do you want your output files in FITS format?
options.FAMIAS_output=false; %Do you want your regularised LSD profiles output ready to load for FAMIAS? Requires overwrite_post_reduction=true and calculate_lsd_profile=true.

%MEGARA CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[final_data,options.objectname,star_database]=megara_control_file(options);
