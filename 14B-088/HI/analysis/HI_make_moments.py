
from spectral_cube import SpectralCube
import numpy as np
from astropy.io import fits

from analysis.paths import fourteenB_HI_data_path
from analysis.constants import (cube_name, mask_name, moment0_name,
                                moment1_name, lwidth_name, skew_name,
                                kurt_name)

'''
Make the first three moments with the pbcov masked cube.

The higher moments (skewness & kurtosis) need a more aggressive mask, made by
make_signal_mask.py, and computed in higher_moments.py.

'''

cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

# Load in source mask
source_mask = fits.getdata(fourteenB_HI_data_path(mask_name))
source_mask = source_mask.astype(np.bool)

cube = cube.with_mask(source_mask)

# Now create the moment 1 and save it. Make a linewidth one too.

moment0 = cube.moment0()
moment0.write(fourteenB_HI_data_path(moment0_name, no_check=True),
              overwrite=True)

moment1 = cube.moment1().astype(np.float32)
moment1.header["BITPIX"] = -32
moment1.write(fourteenB_HI_data_path(moment1_name, no_check=True),
              overwrite=True)

linewidth = cube.linewidth_sigma()
linewidth.write(fourteenB_HI_data_path(lwidth_name, no_check=True),
                overwrite=True)

# Skewness
mom3 = cube.moment(order=3, axis=0)

# Normalize third moment by the linewidth to get the skewness
skew = mom3 / linewidth ** 3
skew.write(fourteenB_HI_data_path(skew_name, no_check=True),
           overwrite=True)

# Kurtosis: Uncorrected
mom4 = cube.moment(order=4, axis=0)
# mom4_resc = (mom4.value - 3) ** 0.25 / 1000.

# Normalize third moment by the linewidth to get the skewness
kurt = mom4 / linewidth ** 4

kurt.write(fourteenB_HI_data_path(kurt_name, no_check=True),
           overwrite=True)
