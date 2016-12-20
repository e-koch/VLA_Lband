
'''
Given a velocity surface, shift all spectra in the cube.
'''

from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
from multiprocessing import Pool

from spectra_shifter import cube_shifter
from paths import fourteenB_HI_data_path
from constants import cube_name, centroidsub_cube_name, moment1_name
from galaxy_params import gal


cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

# Shift with the first moment
mom1 = fits.open(fourteenB_HI_data_path(moment1_name))[0].data * u.m / u.s

pool = Pool(6)

cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_path(centroidsub_cube_name,
                                              no_check=True),
             return_spectra=False, verbose=True, pool=pool)
