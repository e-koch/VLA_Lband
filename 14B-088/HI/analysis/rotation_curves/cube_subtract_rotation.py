
from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
import os
from astropy.utils.console import ProgressBar

'''
Subtract a rotation model from a cube.
'''

# Load in my huge FITS creator
execfile(os.path.expanduser("~/Dropbox/code_development/ewky_scripts/write_huge_fits.py"))


def find_nearest(array, value):
    '''
    http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    '''
    idx = (np.abs(array - value)).argmin()
    return idx

# Set vsys. Using the fit value from DISKFIT
vsys = -180610 * u.m / u.s

data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"

cube = SpectralCube.read(os.path.join(data_path,
                         "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))

# Where's the center?
center_pixel = find_nearest(cube.spectral_axis, vsys)
# In this case, the remaining difference is a minuscule 3 m/s.

model = fits.open(os.path.join(data_path,
                  "diskfit_noasymm_nowarp_output/rad.mod.fits"))

# Now calculate the spectral shifts needed for each pixel
# Assuming that the array shapes for the same (which they are here)
shifts = np.zeros(model[0].data.shape)

posns = np.where(np.isfinite(model[0].data))

# Adjust the header
new_header = cube.header.copy()
# There's a 1 pixel offset
new_header["CRPIX3"] = center_pixel + 1
new_header["CRVAL3"] = (cube.spectral_axis[center_pixel] - vsys).value

# Create the FITS file so we can write 1 spectrum in at a time
new_fitsname = os.path.join(data_path,
                            "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits")

create_huge_fits(cube.shape, new_fitsname, header=new_header)

new_fits = fits.open(new_fitsname, mode='update')

write_every = 1000

for num, (i, j) in enumerate(ProgressBar(zip(*posns))):
    shift = find_nearest(cube.spectral_axis,
                         model[0].data[i, j] * u.m / u.s) - center_pixel
    new_fits[0].data[:, i, j] = np.roll(cube.filled_data[:, i, j].astype(np.float32), shift)

    if num % write_every == 0:
        new_fits.flush()

new_fits.flush()
new_fits.close()
