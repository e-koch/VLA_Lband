
# First Data Release from the M33 L-band VLA Project
# Data presented in from Koch et al. (2018)
# Paper DOI: 
# Code link: https://github.com/ekoch/VLA_Lband/

This is the first data release of the M33 L-band VLA project. This release includes most of the spectral-line data taken in VLA project 14B-088.

Further requests and questions regarding the data can be sent to koch.eric.w@gmail.com.

# HI Data

* `M33_14B-088_HI.clean.image.GBT_feathered.2.6kms.fits` - HI cube with a velocity resolution of 2.6 km/s. **Note:** This is not the cube version presented in the paper linked above. The 200 m/s cube will be released in the next data release.
* `M33_14B-088_pbcov.fits` - The primary beam coverage of the HI data.
* Moment arrays - These are derived from the 200 m/s cube with a primary beam cut-off of 0.5 applied.
* `rotation_curves` - Contains the DISKFIT outputs for fitting a circular rotation curve to the peak velocity surface. `rad.out.csv` contains the circular velocities and errors in each bin. `rad.out.params.csv` contains the fitted disk parameters. FITS files with `fit` in the name are the fitted model or residual compared to fitting a Brandt rotation curve, as described in the paper. `rad.mod.fits` and `rad.res.fits` are the model and residual surfaces directly from the DISKFIT model.
* `stacked_spectra` - A series of stacked spectra over the entire disk, in radial bins, in peak temperature bins, and their respective fits. Note that the stacked profiles are given only for the feathered VLA+GBT data, but the fit tables have fits to both the feathered data and the VLA-only data. Since we consider the former to be more scientifically useful, the VLA-only stacked profiles are not given here. Contact the authors for the VLA-only stacked profiles.
    - `NAME_stacked.fits` - Stacked spectra over the entire disk (complete out to a radius of ~8 kc). There are files for aligning the spectra based on the rotation, centroid and peak velocity surfaces.
    - `hi_gaussian_totalprof_hwhm_fits.csv` - Fit parameters to the previous stacked profiles. See the linked paper for a description of the fit parameters. Includes fits to the VLA-only and VLA+GBT stacked profiles. The latter have "feath" in their column names.
    - `NAME_100pc.fits` - Stacked spectra in 100 pc radial bins (in the frame of the galaxy). Files includes stacks over the entire radial bin and those limited to the northern and southern halves, respectively.
    - `hi_gaussian_hwhm_totalprof_fits_radial.csv` - Fits to the radially-stacked profiles. See the linked paper for a description of the fit parameters. Includes fits to the VLA-only and VLA+GBT stacked profiles. The latter have "feath" in their column names.
    - `NAME_5percentile.fits` - Stacks in 5-percentile intervals of the peak temperature distribution.
    -  `hi_hwhm_totalprof_fits_peak_5percentile_feather.csv` - Fits to the peak temperature stacked profiles. See the linked paper for a description of the fit parameters. Only includes fits to the VLA+GBT stacked profiles.

# OH Data

Note that data for the separately-published OH-1665 maser detection can be found [here](http://apps.canfar.net/storage/list/ekoch/M33/VLA/Koch2018_OH_maser).

Naturally-weighted data cubes with 1.5 km/s velocity channels are given for each of the 4 transitions.  These are dirty images with no deconvolution applied. There is only one clear detection in the 1665 MHz line (see above).

# RRL Data

Naturally-weighted data cubes with 10 km/s velocity channels are given for 5 Hydrogen radio-recombination lines. These are dirty images with no deconvolution applied. There are no clear detections.

Note that some of the cubes show artifacts from unflagged RFI in some pointings. These are left as is since the field with the highest coverage, which covers the giant HII region NGC 604, was found to show no emission, making it unlikely other regions would since it is one of the brightest Halpha sources in the galaxy.

* `HRL_14B-088_combined.fits` - A co-added version of the individual lines smoothed to match the lowest resolution line (H152alp).
* `NAME_14B-088_dirty.image.fits` - Primary-beam, naturally-weighted cube for each line.
* `NAME_14B-088_dirty.image.smooth.fits` - Version smoothed to a common beam for use in the co-added version.
