# MEGARA 2.2

MEGARA is a MATLAB pipeline for reducing and post-processing HERCULES echelle spectra taken on the CCD486.

## Setup

1. Remove older MEGARA versions and their subfolders from the MATLAB path, then add the complete `Version_2_2` tree.
2. Place `pointer_raw_data_directory.m` in the raw-data root directory.
3. Place `pointer_reduced_data_directory.m` in the reduced-data root directory.
4. Add both data root directories to the MATLAB path.

MEGARA uses these pointer files to locate its input and output directories.

Keep raw observations and reduced products outside the `Version_2_2` folder.
Copy the bundled tutorial data to the raw-data root before running it.

## Running MEGARA

Edit and run:

```text
User_codes/Run_megara.m
```

Set `options.objectname` and review the options in that file before running. For an existing reduction, set `options.skip_reduction=true` to run post-reduction only.

MEGARA reuses its MATLAB process pool between observing folders. A two-folder
CCD486 HD110379 check found identical February extracted and processing files
with a reused pool and a fresh pool. That check used cached flat and ThAr
calibrations; a calibration-rebuild transition has not been tested yet.

## Main outputs

Outputs are stored in `Reduced_Data/<target>/`:

- `reduced frames/J*.mat` — extracted order spectra and observation metadata.
- `reduced frames/processing_files/J*_prc.mat` — extraction products.
- `final_data_<target>.mat` — merged, post-reduced spectra.
- `regularised_lsd_profile_<target>.mat` — regularised LSD profiles.
- `supplementary_data/diagnostics_<target>.mat` — continuum-fitting and
  RV/vsini diagnostic structures.
- `supplementary_data/` — databases, cleanup reports, cached continuum fits
  and other supporting information.
- Optional ASCII, FITS, and FAMIAS products.

New final-data files contain a variable named `final_data`. Readers also accept
the older `finaldata` variable. An existing file is left intact if writing its
replacement fails.

The wavelength columns in `J*.mat` are not barycentrically shifted. The correction is stored as `bcorr` and is applied during post-reduction when `options.apply_barycentric_correction=true`.

The barycentric correction is calculated entirely in MATLAB using the
bundled stellar catalogue, DE405 ephemeris and HERCULES observatory model.
No external program or C compilation is required. If a target name is not
identified, MEGARA warns and uses the coordinates in the FITS header without
parallax, proper motion or radial velocity.

## Code layout

- `User_codes/` — run scripts and target databases.
- `Pointer_files/` — templates for the raw- and reduced-data pointer files.
- `Reduction_code/` — main pipeline.
- `Other_useful_code/` — supplementary tools and diagnostics.

The bundled `MEGARA 2.2 Manual.pdf` retains the MEGARA 2.1 technical content;
its cover version and abstract have been updated to 2.2. The installation
instructions in this package describe 2.2-specific setup.

## Stellar extraction

New reductions register each stellar frame's extraction to its monthly master
flat. This requires re-reducing the raw stellar frames; post-reduction alone
does not change an existing extraction. The per-frame `trace_registration`
diagnostic is saved in each reduced frame.

