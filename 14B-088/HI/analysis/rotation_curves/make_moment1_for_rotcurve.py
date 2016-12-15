
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
from lmfit import minimize, Parameters

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

pool.close()
pool.join()

for out in output:
    peakvels[out[1], out[2]] = out[0]

peakvels[peakvels == 0.0 * u.m / u.s] = np.NaN * u.m / u.s

peak_intens_name = "M33_14B-088_HI.clean.image.pbcov_gt_{}.ellip_mask.peakvels.fits".formt(pb_lim)
peakvels = peakvels.astype(np.float32)
# peakvels = spectral_peakintensity(subcube).astype(np.float32)
peakvels.header["BITPIX"] = -32
peakvels.write(fourteenB_HI_data_path(peak_intens_name, no_check=True),
               overwrite=True)


# Trying Gauss-Hermite fits. The basic outline of the lmfit code is from:
# http://research.endlessfernweh.com/curve-fitting/

gauss_herm_center = Projection(np.zeros(subcube.shape[1:]),
                               wcs=cube.wcs.celestial,
                               unit=cube.spectral_axis.unit)


def gaussfunc_gh(paramsin, x):

    try:
        amp = paramsin['amp'].value
        center = paramsin['center'].value
        sig = paramsin['sig'].value
        skew = paramsin['skew'].value
        kurt = paramsin['kurt'].value
    except:
        amp, center, sig, skew, kurt = paramsin

    c1 = -np.sqrt(3)
    c2 = -np.sqrt(6)
    c3 = 2 / np.sqrt(3)
    c4 = np.sqrt(6) / 3.
    c5 = np.sqrt(6) / 4.

    term1 = (x - center) / sig

    gaustot_gh = amp * np.exp(-.5 * term1 ** 2) * \
        (1 + skew * (c1 * term1 + c3 * term1 ** 3) +
         kurt * (c5 + c2 * term1**2 + c4 * term1**4))
    return gaustot_gh


def herm_gauss_fitting(posn):
    y, x = posn
    spec = subcube[:, y, x].value
    spec_axis = subcube.spectral_axis.value
    chanwidth = np.abs(spec_axis[1] - spec_axis[0])

    p_gh = Parameters()
    p_gh.add('amp', value=spec.max(), vary=True)
    p_gh.add('center', value=spec_axis[spec.argmax()], min=np.min(spec_axis),
             max=np.max(spec_axis))
    p_gh.add('sig', value=30 * chanwidth, min=chanwidth, max=None)
    p_gh.add('skew', value=0, vary=True, min=None, max=None)
    p_gh.add('kurt', value=0, vary=True, min=None, max=None)

    def gausserr_gh(p, x, y):
        return gaussfunc_gh(p, x) - y

    fitout_gh = minimize(gausserr_gh, p_gh, args=(spec_axis, spec))

    fit_gh = gaussfunc_gh(fitout_gh.params, spec_axis)

    verbose = False
    if verbose:
        import matplotlib.pyplot as p
        p.plot(spec_axis, fit_gh)

    return spec_axis[fit_gh.argmax()], y, x


pool = Pool(6)
output_hg = pool.map(herm_gauss_fitting, izip(*posns))

pool.close()
pool.join()

for out in output_hg:
    gauss_herm_center[out[1], out[2]] = out[0] * u.m / u.s

gauss_herm_center[gauss_herm_center == 0.0 * u.m / u.s] = np.NaN * u.m / u.s

gauss_herm_name = "M33_14B-088_HI.clean.image.pbcov_gt_{}.ellip_mask.gh_peak.fits".format(pb_lim)
gauss_herm_center = gauss_herm_center.astype(np.float32)
gauss_herm_center.header["BITPIX"] = -32
gauss_herm_center.write(fourteenB_HI_data_path(gauss_herm_name, no_check=True),
                        overwrite=True)
