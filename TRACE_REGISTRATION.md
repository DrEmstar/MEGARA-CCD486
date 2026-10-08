# Stellar extraction registration in MEGARA 2.1

Registration requires a preliminary extraction and therefore adds processing time. New stellar reductions register the extraction aperture to the science frame before applying the monthly flat's spatial weights. This addresses continuum ripples caused by relative motion between the stellar illumination and the monthly master trace. The monthly flat, its response, its profile weights and the ThAr calibration are not modified.

## Algorithm

1. Extract the star and fitted background on the original monthly trace. Apply the existing cosmic-ray filter to the star.
2. Measure positive-flux spatial centroids of the background-subtracted star and the background-subtracted master-flat profiles, in the same flat-selected aperture and at identical detector columns.
3. Reject nonfinite, poorly illuminated and outlying columns. Require at least 100 usable columns per order. Combine reliable per-order offsets into one robust detector-wide shift; orders with discrepant shifts are excluded.
4. Re-extract the original star and fitted background using the shifted trace, and apply the existing cosmic-ray filter. Continue with the original flat weighting, flat-fielding, normalisation, uncertainties and merge weights.

The estimator uses all available orders, rather than only the five blue orders used in the initial diagnostic. It requires at least three accepted orders, rejects within-order centroid scatter above 1 pixel, inter-order robust scatter above 0.5 pixel, and offsets exceeding 5 pixels. These are conservative engineering bounds, not calibrated probability thresholds. Failed registration retains the original trace and emits a warning. If shifting an individual order changes its extracted detector columns, that order retains its original extraction, preserving wavelength and flat indexing.

Both the reduced exposure and its processing MAT file contain `trace_registration`: method version, status, measured/applied shift, per-order offsets, valid-column counts, accepted orders, scatter and any per-order fallbacks. Order indices in these diagnostics refer to the monthly trace numbering, before conversion to absolute echelle orders. Offsets are science minus monthly-flat centroids in detector pixels. A positive shift increases the extraction trace's detector row.

## Configuration and existing reductions

Registration is on by default. For a controlled comparison, set `stellar_trace_registration_enabled=false` in `Reduction_code/reduce_hercules_frames.m`. The runtime flag is set after saving/loading the monthly flat cache, so cached calibrations do not determine this choice. Direct callers may set the same field on `allorders`.

Existing reduced files are still skipped by the pipeline. Updating the code or rerunning post-reduction alone does not apply this correction to them. To re-reduce an existing target, first archive its old extracted products and ensure they are no longer in either the searched raw folder or reduced target folder; preserve the raw FITS and monthly calibration files. Then rerun raw reduction and downstream continuum processing. Do not mix fixed-trace and registered extractions in a comparison without checking the saved metadata.

## Evidence and limits

In the original 20-exposure HD172555 test (orders 142–146, 3900–4000 Å), the median fitted ripple semi-amplitude relative to low-ripple controls fell from 1.55% to 0.23%; all 12 exposures initially above 1% improved. The strongest two fell from 6.50% and 6.36% to 0.24% and 0.38%. Those numbers describe the diagnostic five-order estimator, not an absolute continuum calibration or a guarantee for every target.

Six initially weak-ripple exposures did not improve, including a 2018 exposure. Follow-up showed sensitivity to spatial profile weighting and reference choice. Registration alone does not correct all illumination-profile changes. The wavelength solution remains unchanged: apparent Ca II profile shifts of up to approximately 0.18 km/s depend on fitting window and continuum model. Precision radial-velocity and line-variability analyses require additional calibration/profile validation. No claim of full-spectrum continuum or velocity accuracy is made.

Run the deterministic regression checks in MATLAB with:

```matlab
addpath('tests');
test_stellar_trace_registration
```

The tests cover positive/negative offsets, flat-centroid bias, order outliers, unavailable signal/calibration, the off switch, excessive offsets, nonfinite values, and detector-edge fallback.

## Implementation validation (October 2026)

The production estimator was checked against the same 20 raw exposures and their respective monthly flats. Four frames (J8303028, J0845038, J0862142 and J0867013) exercised full-resolution extraction across all detector orders. For the other 16, validation sampled every sixteenth detector column for the registration estimate, then re-extracted orders 142–146 at full resolution. Production always uses all available columns. The four full-resolution runs required no detector-edge fallbacks.

Using identical downstream continuum/merge processing, the median differential ripple semi-amplitude was 1.547% before and 0.243% after. All 12 initially above 1% improved; 14 of 20 improved overall. J0867013 changed from 6.503% to 0.247%; J0867014 from 6.362% to 0.354%. The median high-frequency noise-proxy ratio was 1.004. Weak cases can worsen: J8303028 changed from 0.233% to 0.722%, and J0760056 from 0.304% to 1.025%. A detector-wide translation can differ from a local blue-order offset; these checks do not establish accuracy throughout the spectrum.

Deterministic regression tests and MATLAB code-analysis checks passed. The disabled path reproduces the previous extracted arrays exactly in the regression test. The intentional detector-edge failure test confirms retention of the original order without changing its detector-column indexing.
