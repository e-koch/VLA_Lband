
from spectral_cube import SpectralCube
import os
import pyregion
import astropy.units as u
from astropy.io import fits


path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging"

cube_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"

cube = SpectralCube.read(os.path.join(path, cube_name))

# Mask at 3 sigma
cube = cube.with_mask(cube > 3 * 1.8 * u.mJy)

# Now cut to the elliptical region to remove all bkg regions
region = pyregion.open(os.path.expanduser("~/Dropbox/code_development/VLA_Lband"
                                          "/14B-088/HI/analysis/mom1_rotcurve_mask.reg"))
subcube = cube.subcube_from_ds9region(region)

# Now create the moment 1 and save it. Make a linewidth one too.
# DISKFIT has issues with float64, so convert to float32 then save

moment0_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3.ellip_mask.mom0.fits"
moment0 = subcube.moment0()
moment0.write(os.path.join(path, moment0_name))

moment1_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3.ellip_mask.mom1.fits"
moment1 = subcube.moment1().astype(np.float32)
moment1.header["BITPIX"] = -32
moment1.write(os.path.join(path, moment1_name))

lwidth_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3.ellip_mask.linewidth.fits"
linewidth = subcube.linewidth_sigma()
linewidth.write(os.path.join(path, lwidth_name))
