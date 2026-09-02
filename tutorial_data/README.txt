MEGARA 2.1 tutorial data

This folder contains the raw FITS frames for the HD 48501 tutorial described
in Appendix A of the MEGARA 2.1 manual.

Before running the tutorial:

1. Copy the complete tutorial_data folder to your raw-data root.
2. Keep the raw-data and reduced-data roots outside the MEGARA program folder.
3. Set options.reduction_folders = {'tutorial_data'} in User_codes/Run_megara.m.

Do not run the reduction in this release copy. MEGARA writes calibration and
classification records into the observing-run directory, so work from the
copy in your raw-data root instead.

Only original raw FITS frames are supplied. Cached calibration, classification
and reduction products from earlier MEGARA versions have intentionally been
excluded so that the tutorial exercises the MEGARA 2.1 workflow.
