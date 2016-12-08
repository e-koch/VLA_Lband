
from astropy.io import fits
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
import os
from astropy.utils.console import ProgressBar
import pyregion
from astropy.table import Table
from turbustat.statistics.stats_utils import fourier_shift

# Import from above.
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)
from rotation_curves.vrot_fit import return_smooth_model
from paths import fourteenB_HI_data_path, iram_co21_data_path
from galaxy_params import gal

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
table_name = fourteenB_HI_data_path("diskfit_noasymm_noradial_nowarp_output/rad.out.csv")
tab = Table.read(table_name)

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
# Somehow blank keywords are generated when reading the cube in... something
# to do with dropping the Stokes axis is my guess.
del cube._header[""]

# Where's the center?
# center_pixel = find_nearest(cube.spectral_axis, vsys)
# In this case, the remaining difference is a minuscule 3 m/s.

model = return_smooth_model(tab, cube.header, gal)

# Now calculate the spectral shifts needed for each pixel
# Assuming that the array shapes for the same (which they are here)
shifts = np.zeros(model.shape)

posns = np.where(np.isfinite(model))

# Adjust the header
new_header = cube.header.copy()
# There's a 1 pixel offset
# new_header["CRPIX3"] = center_pixel + 1
new_header["CRVAL3"] = new_header["CRVAL3"] - gal.vsys.to(u.m / u.s).value

# Create the FITS file so we can write 1 spectrum in at a time
print("Making the empty FITS file")
new_fitsname = iram_co21_data_path("m33.co21_iram.rotsub.fits", no_check=True)

if os.path.exists(new_fitsname):
    print("Removing old rotsub version.")
    os.system("rm {}".format(new_fitsname))

create_huge_fits(new_fitsname, new_header)

new_fits = fits.open(new_fitsname, mode='update')

write_every = 20000

vel_res = cube._pix_size_slice(0) * u.m / u.s

vsys = gal.vsys.to(u.m / u.s)

print("And here we go!")
for num, (i, j) in enumerate(ProgressBar(zip(*posns))):
    # Don't bother rolling if there's nothing there
    if not np.isfinite(cube.filled_data[:, i, j]).any():
        continue

    shift = ((model[i, j] * u.m / u.s - vsys) / vel_res).value

    new_fits[0].data[:, i, j] = \
        fourier_shift(cube.filled_data[:, i, j], shift)

    if num % write_every == 0:
        new_fits.flush()

new_fits.flush()
new_fits.close()
