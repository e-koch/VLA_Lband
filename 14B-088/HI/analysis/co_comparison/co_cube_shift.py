
'''
Shift the CO data with different velocity surfaces.
'''

from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
from multiprocessing import Pool
from reproject import reproject_interp

from spectra_shifter import cube_shifter
from paths import fourteenB_HI_data_path, iram_co21_data_path
from constants import moment1_name
from galaxy_params import gal


cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
del cube._header[""]

# Shift with the first moment
mom1 = fits.open(fourteenB_HI_data_path(moment1_name))[0]

# Reproject the centroid map
mom1_reproj = reproject_interp(mom1, cube.wcs.celestial,
                               shape_out=cube.shape[1:])[0] * u.m / u.s

pool = Pool(6)

cube_shifter(cube, mom1_reproj, gal.vsys, save_shifted=True,
             save_name=iram_co21_data_path("m33.co21_iram.hi_centroid_corrected.fits",
                                           no_check=True),
             return_spectra=False, verbose=True, pool=pool)

pool.close()
pool.join()
