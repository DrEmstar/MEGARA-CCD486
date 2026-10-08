# MEGARA 2.1

MEGARA is a MATLAB pipeline for reducing and post-processing HERCULES echelle spectra taken on the CCD486.

## Setup

1. Add the complete `Version_2_1` folder tree to the MATLAB path.
2. Place `pointer_raw_data_directory.m` in the raw-data root directory.
3. Place `pointer_reduced_data_directory.m` in the reduced-data root directory.
4. Add both data root directories to the MATLAB path.

MEGARA uses these pointer files to locate its input and output directories.

## Running MEGARA

Edit and run:

```text
User_codes/Run_megara.m
```

Set `options.objectname` and review the options in that file before running. For an existing reduction, set `options.skip_reduction=true` to run post-reduction only.

## Stellar extraction alignment

New reductions align the stellar extraction to the science frame relative to its
monthly master flat. The measured offset and fallback diagnostics are saved with
each exposure. See `TRACE_REGISTRATION.md` for the method, off switch, instructions
for re-reducing existing data, validation and wavelength/profile limitations.

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
- `Deicoon/` — LSD resources.
- `Other_useful_code/` — supplementary tools and diagnostics.

The full `MEGARA 2.1 Manual.pdf` is included in the repository root.

## Release contents

This branch is the CCD486 MEGARA 2.1 release. It is maintained in parallel
with the existing MEGARA 1.6 Mac and Unix branches; those legacy branches are
unchanged. The bundled `tutorial_data` directory contains the HD 48501 files
used by the tutorial in the manual.
