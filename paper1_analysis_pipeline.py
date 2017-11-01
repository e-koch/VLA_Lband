
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

import os
import sys
from subprocess import call

from paths import (c_hi_analysispath, fourteenB_HI_data_path,
                   fourteenB_HI_data_wGBT_path, fourteenB_HI_file_dict,
                   fourteenB_wGBT_HI_file_dict)
from analysis.constants import moment1_name

os.sys.path.insert(0, c_hi_analysispath(""))
# os.chdir(c_hi_analysispath(""))

#################
# Initial masking.
#################

# Mask w/ pbcov limit.
execfile(c_hi_analysispath("pbcov_masking.py"))

# Create the source mask
execfile(c_hi_analysispath("make_signal_mask.py"))

# Create moment arrays
execfile(c_hi_analysispath("HI_make_moments.py"))

# Make some moment figures
execfile(c_hi_analysispath("HI_moment_figures.py"))

#################
# Rotation curves
#################

# Find the rotation curve
# Create a centroid map w/ a radial cutoff
execfile(c_hi_analysispath("rotation_curves/make_moment1_for_rotcurve.py"))

# Now run diskfit.

# First using the centroids
diskfit_script = c_hi_analysispath("rotation_curves/run_diskfit.py")
diskfit_params = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_nowarp_noradial_noasymm.inp")
output_path = fourteenB_HI_data_path("", no_check=True)

call("python {0} {1} {2} {3}".format(diskfit_script, diskfit_params,
                                     output_path,
                                     fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.5.ellip_mask.mom1.fits")),
     shell=True)

# And then using peak velocity
diskfit_params = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_peakvels_nowarp_noradial_noasymm.inp")
output_path = fourteenB_HI_data_path("", no_check=True)

call("python {0} {1} {2} {3}".format(diskfit_script, diskfit_params,
                                     output_path,
                                     fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.5.ellip_mask.mom1.fits")),
     shell=True)

# And with the peak from a Gauss-Hermite polynomial
diskfit_params = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_ghfit_nowarp_noradial_noasymm.inp")
output_path = fourteenB_HI_data_path("", no_check=True)

call("python {0} {1} {2} {3}".format(diskfit_script, diskfit_params,
                                     output_path,
                                     fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.5.ellip_mask.mom1.fits")),
     shell=True)

# Compare the three rotation curves.
execfile(c_hi_analysispath("rotation_curves/rotation_curve_comparison.py"))

# Fit the rotation curve and produce a smooth model
call("python {0} {1}".format(c_hi_analysispath("rotation_curves/vrot_fit.py"),
                             fourteenB_HI_data_path("diskfit_noasymm_noradial_nowarp_output")), shell=True)

# Residual velocity figures
execfile(c_hi_analysispath("rotation_curves/residual_velocity_surface_figure.py"))

# Subtract the rotation curve from the cube and save.
execfile(c_hi_analysispath("rotation_curves/cube_subtract_rotation.py"))

# Rotation subtracted slab figures
execfile(c_hi_analysispath("rotation_curves/rotsub_channels_figure.py"))

# Movie! Turned off by default.
# execfile(c_hi_analysispath("rotsub_channels_movie.py"))

# Save the new galactic parameters value as a table for the paper
call("python {0} {1}".format(c_hi_analysispath("rotation_curves/save_gal_params_table.py"),
                             fourteenB_HI_data_path("diskfit_noasymm_noradial_nowarp_output")), shell=True)

#################
# HI analysis
#################

# Surface Density profiles
execfile(c_hi_analysispath("HI_radial_profile.py"))

# Velocity dispersion profile
execfile(c_hi_analysispath("HI_veldisp_profile.py"))

# Plots of skewness and kurtosis
execfile(c_hi_analysispath("higher_moments.py"))

# Total galaxy profile.
execfile(c_hi_analysispath("total_profiles.py"))

# Run filament analysis on zeroth moment
execfile(c_hi_analysispath("run_filfinder.py"))

########
# Other
########

# Make a UV plane plot
execfile(c_hi_analysispath("uv_plots/channel_1000_uvplot.py"))
