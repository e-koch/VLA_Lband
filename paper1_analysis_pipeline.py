
'''
Run entire analysis. The cube should already have been exported into a FITS
file from CASA.

This script includes deriving the signal mask and moment arrays for the VLA
and VLA + GBT HI cubes.  The feathering and comparisons should already have
been run.

Paper 1 investigates the ensemble of HI properties:
 - using diskfit to fit rotation models to the VLA and VLA + GBT peak velocity
 - relations to the skewness and kurtosis
 - effect of centering with the rotation model, centroid, and peak velocity
 - line shapes of stacked profiles
 - comparing estimates of the line width
'''

from paths import c_hi_analysispath

#################
# Initial masking.
#################

execfile(c_hi_analysispath("cube_pipeline.py"))

# Make some moment figures
execfile(c_hi_analysispath("HI_moment_figures.py"))

#################
# Rotation curves
#################

# There needs to be a manual step here where the centroid, rotation and peak
# velocity maps are re-saved in a "moments_for_diskfit" as 32-bit floats.

execfile(c_hi_analysispath("rotation_curves/fit_rotation_curves.py"))
execfile(c_hi_analysispath("rotation_curves/fit_rotation_curves_feathered.py"))

# Compare the three rotation curves.
# execfile(c_hi_analysispath("rotation_curves/rotation_curve_comparison.py"))

# Residual velocity figures
execfile(c_hi_analysispath("rotation_curves/residual_velocity_surface_figure.py"))

# Figures
execfile(c_hi_analysispath("rotation_curves/rotation_curve_figures.py"))

# Rotation subtracted slab figures
# execfile(c_hi_analysispath("rotation_curves/rotsub_channels_figure.py"))

# Movie! Turned off by default.
# execfile(c_hi_analysispath("rotsub_channels_movie.py"))

# Make velocity subtracted cubes
execfile(c_hi_analysispath("cube_shift_spectra.py"))
execfile(c_hi_analysispath("cube_shift_spectra_feathered.py"))

#################
# HI analysis
#################

# Data example figures
execfile(c_hi_analysispath("HI_moment_figures.py"))
execfile(c_hi_analysispath("HI_example_spectra_figure.py"))

# Masses
execfile(c_hi_analysispath("HI_masses.py"))

# Surface Density profiles
execfile(c_hi_analysispath("HI_radial_profile.py"))

# Plots of skewness and kurtosis
execfile(c_hi_analysispath("higher_moments.py"))
execfile(c_hi_analysispath("skew_kurt_properties.py"))

# Total galaxy profile.
# XXX Update this path!
execfile(c_hi_analysispath("total_profiles.py"))

# Run filament analysis on zeroth moment
# execfile(c_hi_analysispath("run_filfinder.py"))

# Stacking and analysis
execfile(c_hi_analysispath("HI_radial_stacking.py"))
execfile(c_hi_analysispath("HI_radial_stacking_feathered.py"))
execfile(c_hi_analysispath("HI_stacking_modeling.py"))
execfile(c_hi_analysispath("HI_radial_stacking_NS_peak_modeling.py"))

# Velocity dispersion profile
execfile(c_hi_analysispath("HI_veldisp_profile.py"))

# PV-slices
execfile(c_hi_analysispath("HI_pvslices.py"))
execfile(c_hi_analysispath("HI_pvslices_figures.py"))
execfile(c_hi_analysispath("HI_pvslices_nplume.py"))
execfile(c_hi_analysispath("HI_pvslices_thin.py"))
execfile(c_hi_analysispath("HI_pvslices_thin_figures.py"))


########
# Other
########

# Make a UV plane plot
# execfile(c_hi_analysispath("uv_plots/channel_1000_uvplot.py"))
