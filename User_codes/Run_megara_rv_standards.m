%Reduce the high-vsini validation targets with MEGARA 2.1.

clearvars
clc

%Use the archive CCD data for MEGARA 2.1, not the new CCD directory.
archive_raw_data_directory=['/Volumes/TOSHIBA EXT/Astronomy_Data/' ...
    'Mt John Data/Jfiles_bydate'];
if ~exist(fullfile(archive_raw_data_directory, ...
        'pointer_raw_data_directory.m'),'file')
    error('The MEGARA 2.1 archive is not mounted at %s.', ...
        archive_raw_data_directory)
end
original_matlab_path=path;
path_cleanup=onCleanup(@() path(original_matlab_path));
addpath(archive_raw_data_directory,'-begin')
rehash path

%Stellar parameters used to select the synthetic spectrum and LSD mask.
%Update these values here if preferred standard parameters are adopted.
standards={ ...
    'hd_149989', 7250, 4.0, 136; ...
    'hd_224945', 7250, 4.0, 55};

%Shared MEGARA options
options=struct();
options.reduction_folders={};

options.star_census=true;
options.skip_reduction=false;
options.skip_post_reduction=false;
options.skip_manual_file_check=true;

options.chunked_post_reduction=false;
options.post_reduction_chunk_size=500;
options.extended_wavelength_range=false;

options.blue_data_chop_value=300;
options.apply_barycentric_correction=true;

options.apply_continuum_fit=true;
options.manual_continuum_fit=false;
options.template_continuum_fit=false;
options.automatic_continuum_fit=true;
options.order_merge=true;

options.overwrite_full_data=true;
options.overwrite_post_reduction=true;
options.calculate_lsd_profile=true;
options.radial_velocity_measurement=true;
options.extend_velocities=false;

options.suppress_figures=true;
options.ascii_output=false;
options.fits_output=false;
options.FAMIAS_output=false;
options.order_rv_diagnostic=true;

%Run each standard separately because the stellar parameters differ
standard_results=struct();

for standard_index=1:size(standards,1)

    options.objectname=standards{standard_index,1};
    options.teff=standards{standard_index,2};
    options.logg=standards{standard_index,3};
    options.vsini=standards{standard_index,4};

    fprintf('\nReducing RV validation target %s (%d of %d)\n', ...
        options.objectname,standard_index,size(standards,1))

    %Allow any failure to reach the command window with its full stack.
    [final_data,reduced_objectname,star_database]= ...
        megara_control_file(options);
    standard_results.(options.objectname).final_data=final_data;
    standard_results.(options.objectname).objectname=reduced_objectname;
    standard_results.(options.objectname).star_database=star_database;
    standard_results.(options.objectname).completed=true;
    standard_results.(options.objectname).error='';

end

completed_standards=fieldnames(standard_results);
number_completed=0;
for standard_index=1:numel(completed_standards)
    number_completed=number_completed+ ...
        standard_results.(completed_standards{standard_index}).completed;
end

fprintf('\nCompleted %d of %d RV-validation reductions\n', ...
    number_completed,size(standards,1))

%Create the literature and method comparison from all available results.
summarise_rv_standard_comparison
