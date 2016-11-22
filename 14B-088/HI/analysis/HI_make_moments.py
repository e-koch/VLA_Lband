
from spectral_cube import SpectralCube
import os
import astropy.units as u
import numpy as np

from analysis.paths import fourteenB_HI_data_path

'''
Make the first three moments with the pbcov masked cube.

The higher moments (skewness & kurtosis) need a more aggressive mask, made by
make_signal_mask.py, and computed in higher_moments.py.

'''
cube_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"

cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

# Mask at 3 sigma
cube = cube.with_mask(cube > 3 * 1.8 * u.mJy)

# Now create the moment 1 and save it. Make a linewidth one too.

moment0_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits"
moment0 = cube.moment0()
moment0.write(fourteenB_HI_data_path(moment0_name, no_check=True))

moment1_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom1.fits"
moment1 = cube.moment1().astype(np.float32)
moment1.header["BITPIX"] = -32
moment1.write(fourteenB_HI_data_path(moment1_name, no_check=True))

lwidth_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.lwidth.fits"
linewidth = cube.linewidth_sigma()
linewidth.write(fourteenB_HI_data_path(lwidth_name, no_check=True))
