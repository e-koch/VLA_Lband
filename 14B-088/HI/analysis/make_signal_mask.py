
from astropy.io import fits
import astropy.units as u
from spectral_cube import SpectralCube
from signal_id import RadioMask, Noise
from scipy import ndimage as nd
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.utils.console import ProgressBar
import skimage.morphology as mo
import numpy as np
from spectral_cube.cube_utils import average_beams

from basics.utils import sig_clip

from analysis.paths import fourteenB_HI_data_path

'''
Create a signal mask for the 14B-088 cube
'''


cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))
# cube = SpectralCube.read("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits")

# noise = Noise(cube)

# scale = noise.scale

scale = sig_clip(cube[-1].value, nsig=10)

pixscale = proj_plane_pixel_scales(cube.wcs)[0]

cube = cube.with_mask(cube > 5 * scale * u.Jy)

# Want to smooth the mask edges
mask = cube.mask.include()

kernel = average_beams(cube.beams).as_tophat_kernel(pixscale)

for i in ProgressBar(xrange(mask.shape[0])):
    mask[i] = nd.binary_opening(mask[i], kernel)
    mask[i] = nd.binary_closing(mask[i], kernel)

    mask[i] = nd.binary_dilation(mask[i], mo.disk(10))

    # Add peak brightness in region check to remove spurious small features

new_header = cube.header.copy()
new_header["BUNIT"] = ""
new_header["BITPIX"] = 8

mask_hdu = fits.PrimaryHDU(mask.astype('>i2'), header=new_header)
mask_hdu.writeto(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked_source_mask.fits",
                                        no_check=True),
                 clobber=True)
# mask_hdu.writeto("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub_source_mask.fits",
#                  clobber=True)
