
'''
Run entire analysis. The cube should already have been exported into a FITS
file from CASA.
'''

import os
import sys
from subprocess import call

from analysis.paths import c_hi_analysispath, fourteenB_HI_data_path
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
                                     fourteenB_HI_data_path(moment1_name)),
     shell=True)

# And then using peak velocity
diskfit_params = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_peakvels_nowarp_noradial_noasymm.inp")
output_path = fourteenB_HI_data_path("", no_check=True)

call("python {0} {1} {2} {3}".format(diskfit_script, diskfit_params,
                                     output_path,
                                     fourteenB_HI_data_path(moment1_name)),
     shell=True)


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

#################
# CO comparison
#################

# Reproject the integrated CO intensity onto the HI grid.
execfile(c_hi_analysispath("co_comparison/co_reproject.py"))

# Subtract the rotation curve from the CO data.
execfile(c_hi_analysispath("co_comparison/co_subtract_rotation.py"))

# Plot where cprops thinks there are clouds.
execfile(c_hi_analysispath("co_comparison/cloud_catalog.py"))

# Radial profiles w/ CO
execfile(c_hi_analysispath("co_comparison/co_radial_profile.py"))

execfile(c_hi_analysispath("co_comparison/co_veldisp_profile.py"))

# Adaptive thresholding mask comparison
execfile(c_hi_analysispath("co_comparison/co_vs_hi_boundaries.py"))

# Rotation subtracted total spectra
execfile(c_hi_analysispath("co_comparison/total_rotsub_profiles.py"))

# HI and CO spectra at the GMC locations
execfile(c_hi_analysispath("co_comparison/co_hi_cloud_spectra.py"))
