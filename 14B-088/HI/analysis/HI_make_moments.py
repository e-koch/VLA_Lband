
from spectral_cube import SpectralCube
import numpy as np
from astropy.io import fits

from analysis.paths import fourteenB_HI_data_path

'''
Make the first three moments with the pbcov masked cube.

The higher moments (skewness & kurtosis) need a more aggressive mask, made by
make_signal_mask.py, and computed in higher_moments.py.

'''
cube_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"

cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

# Load in source mask
source_mask = fits.getdata(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked_source_mask.fits"))
source_mask = source_mask.astype(np.bool)

cube = cube.with_mask(source_mask)

# Now create the moment 1 and save it. Make a linewidth one too.

moment0_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits"
moment0 = cube.moment0()
moment0.write(fourteenB_HI_data_path(moment0_name, no_check=True),
              overwrite=True)

moment1_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom1.fits"
moment1 = cube.moment1().astype(np.float32)
moment1.header["BITPIX"] = -32
moment1.write(fourteenB_HI_data_path(moment1_name, no_check=True),
              overwrite=True)

lwidth_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.lwidth.fits"
linewidth = cube.linewidth_sigma()
linewidth.write(fourteenB_HI_data_path(lwidth_name, no_check=True),
                overwrite=True)

# Skewness
mom3 = cube.moment(order=3)

# Normalize third moment by the linewidth to get the skewness
skew = mom3 / linewidth ** 3
skew.write(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.skewness.fits", no_check=True),
           overwrite=True)

# Kurtosis: Uncorrected
mom4 = cube.moment(order=4)
# mom4_resc = (mom4.value - 3) ** 0.25 / 1000.

# Normalize third moment by the linewidth to get the skewness
kurt = mom4 / linewidth ** 4

kurt.write(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.kurtosis.fits", no_check=True),
           overwrite=True)
