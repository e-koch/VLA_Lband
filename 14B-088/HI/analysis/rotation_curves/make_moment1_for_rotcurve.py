
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import os
import pyregion
import astropy.units as u
import numpy as np
from astropy.io import fits
from scipy.signal import medfilt
from multiprocessing import Pool
from itertools import izip


from analysis.paths import fourteenB_HI_data_path
from analysis.constants import cube_name, mask_name, pb_lim
from analysis.galaxy_params import gal

cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))
mask = fits.open(fourteenB_HI_data_path(mask_name))[0]

# Apply the source mask
cube = cube.with_mask((mask.data > 0))

# Now cut to the elliptical region to remove all bkg regions
# region = pyregion.open(c_hi_analysispath("rotation_curves/mom1_rotcurve_mask.reg"))
# subcube = cube.subcube_from_ds9region(region)

# Since the parameters of M33 are already well constrained, use a radius
# cut-off based on previous values. The fitting is done out to 10 kpc. A
# small bit is added here to account for small changes in the fit parameters.
radius = gal.radius(header=cube.header)
max_radius = 10.25 * u.kpc
subcube = cube.with_mask(radius < max_radius).minimal_subcube()

# Now create the moment 1 and save it. Make a linewidth one too.
# DISKFIT has issues with float64, so convert to float32 then save

moment0_name = "M33_14B-088_HI.clean.image.pbcov_gt_{}.ellip_mask.mom0.fits".format(pb_lim)
moment0 = subcube.moment0()
moment0.write(fourteenB_HI_data_path(moment0_name, no_check=True),
              overwrite=True)

moment1_name = "M33_14B-088_HI.clean.image.pbcov_gt_{}.ellip_mask.mom1.fits".format(pb_lim)
moment1 = subcube.moment1().astype(np.float32)
moment1.header["BITPIX"] = -32
moment1.write(fourteenB_HI_data_path(moment1_name, no_check=True),
              overwrite=True)

lwidth_name = "M33_14B-088_HI.clean.image.pbcov_gt_{}.ellip_mask.linewidth.fits".format(pb_lim)
linewidth = subcube.linewidth_sigma()
linewidth.write(fourteenB_HI_data_path(lwidth_name, no_check=True),
                overwrite=True)

# Make an array of the velocities of the peak intensities.


# def spectral_peakintensity(cube):
#     """
#     Compute the spectral position of the peak intensities.
#     """

#     def peak_velocity(arr, axis=None, smooth_size=31):
#         argmax = np.argmax(medfilt(arr, smooth_size))
#         return cube.spectral_axis[argmax]

#     return cube.apply_function(peak_velocity, axis=0, projection=True,
#                                unit=cube.spectral_axis.unit)

def peak_velocity(inputs):

    y, x = inputs
    smooth_size = 31
    argmax = np.argmax(medfilt(subcube[:, y, x].value, smooth_size))
    return subcube.spectral_axis[argmax], y, x


peakvels = Projection(np.zeros(subcube.shape[1:]),
                      wcs=cube.wcs.celestial,
                      unit=cube.spectral_axis.unit)

posns = np.where(subcube.mask.include().sum(0) > 0)

pool = Pool(6)
output = pool.map(peak_velocity, izip(*posns))

for out in output:
    peakvels[out[1], out[2]] = out[0]

peakvels[peakvels == 0.0 * u.m / u.s] = np.NaN * u.m / u.s

peak_intens_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.5.ellip_mask.peakvels.fits"
peakvels = peakvels.astype(np.float32)
# peakvels = spectral_peakintensity(subcube).astype(np.float32)
peakvels.header["BITPIX"] = -32
peakvels.write(fourteenB_HI_data_path(peak_intens_name, no_check=True),
               overwrite=True)
