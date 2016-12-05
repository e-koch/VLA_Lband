
from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
import os
from astropy.utils.console import ProgressBar
import pyregion
from turbustat.statistics.stats_utils import fourier_shift

# Import from above.
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
from rotation_curves.vrot_fit import return_smooth_model

'''
DUPLICATED CODE WITH ../cube_subtract_rotation.py! A generalized script is
needed!

Subtract a rotation model from the CO cube.
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

data_path = "/media/eric/Data_3/M33/co21/"

cube = SpectralCube.read(os.path.join(data_path, "m33.co21_iram.fits"))

# Where's the center?
center_pixel = find_nearest(cube.spectral_axis, vsys)
# In this case, the remaining difference is a minuscule 3 m/s.

model = return_smooth_model(cube.header)

# Now calculate the spectral shifts needed for each pixel
# Assuming that the array shapes for the same (which they are here)
shifts = np.zeros(model.shape)

posns = np.where(np.isfinite(model))

# Adjust the header
new_header = cube.header.copy()
# There's a 1 pixel offset
new_header["CRPIX3"] = center_pixel + 1
new_header["CRVAL3"] = (cube.spectral_axis[center_pixel] - vsys).value

# Create the FITS file so we can write 1 spectrum in at a time
print("Making the empty FITS file")
new_fitsname = os.path.join(data_path,
                            "m33.co21_iram.rotsub.fits")

create_huge_fits(new_fitsname, new_header)

new_fits = fits.open(new_fitsname, mode='update')

write_every = 20000

vel_res = cube._pix_size_slice(0) * u.m / u.s

print("And here we go!")
for num, (i, j) in enumerate(ProgressBar(zip(*posns))):
    # Don't bother rolling if there's nothing there
    if not np.isfinite(cube.filled_data[:, i, j]).any():
        continue

    shift = ((model[0].data[i, j] * u.m / u.s - vsys) / vel_res).value

    new_fits[0].data[:, i, j] = \
        fourier_shift(cube.filled_data[:, i, j], shift)

    if num % write_every == 0:
        new_fits.flush()

new_fits.flush()
new_fits.close()
