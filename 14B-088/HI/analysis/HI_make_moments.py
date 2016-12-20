
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import numpy as np
import astropy.units as u
from astropy.io import fits
from scipy.signal import medfilt
from itertools import izip
from multiprocessing import Pool

from analysis.paths import fourteenB_HI_data_path
from analysis.constants import (cube_name, mask_name, moment0_name,
                                moment1_name, lwidth_name, skew_name,
                                kurt_name, peakvels_name)

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


# Peak velocity map
def peak_velocity(inputs):

    y, x = inputs
    smooth_size = 31
    argmax = np.argmax(medfilt(cube[:, y, x].value, smooth_size))
    return cube.spectral_axis[argmax], y, x


peakvels = Projection(np.zeros(cube.shape[1:]),
                      wcs=cube.wcs.celestial,
                      unit=cube.spectral_axis.unit)


posns = np.where(cube.mask.include().sum(0) > 0)

pool = Pool(6)
output = pool.map(peak_velocity, izip(*posns))

pool.close()
pool.join()

for out in output:
    peakvels[out[1], out[2]] = out[0]

peakvels[peakvels == 0.0 * u.m / u.s] = np.NaN * u.m / u.s

peakvels = peakvels.astype(np.float32)
# peakvels = spectral_peakintensity(subcube).astype(np.float32)
peakvels.header["BITPIX"] = -32
peakvels.write(fourteenB_HI_data_path(peakvels_name, no_check=True),
               overwrite=True)
