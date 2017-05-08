
'''
Given a velocity surface, shift all spectra in the cube.

Makes 1 version shifted wrt. the centroids.
Makes another shifted wrt. the peak velocity.

'''

from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
from multiprocessing import Pool

from cube_analysis.spectra_shifter import cube_shifter

from paths import fourteenB_HI_data_path
from constants import (cube_name, centroidsub_cube_name, moment1_name,
                       mask_name, centroidsub_mask_name, peakvels_name,
                       peakvelssub_mask_name, peakvelsub_cube_name)
from galaxy_params import gal


cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

mom1 = fits.open(fourteenB_HI_data_path(moment1_name))[0].data * u.m / u.s
peakvels = fits.open(fourteenB_HI_data_path(peakvels_name))[0].data * u.m / u.s

pool = Pool(6)

# Shift with the first moment
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_path(centroidsub_cube_name,
                                              no_check=True),
             return_spectra=False, verbose=True, pool=pool)

# Also shift the signal mask to match those shifted here.

mask_cube = SpectralCube.read(fourteenB_HI_data_path(mask_name))

cube_shifter(mask_cube, mom1, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_path(centroidsub_mask_name,
                                              no_check=True),
             return_spectra=False, verbose=True, pool=pool,
             is_mask=True)


# Shift wrt. peak velocities
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_path(peakvelsub_cube_name,
                                              no_check=True),
             return_spectra=False, verbose=True, pool=pool)

# Also shift the signal mask to match those shifted here.

mask_cube = SpectralCube.read(fourteenB_HI_data_path(mask_name))

cube_shifter(mask_cube, peakvels, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_path(peakvelssub_mask_name,
                                              no_check=True),
             return_spectra=False, verbose=True, pool=pool,
             is_mask=True)

pool.close()
pool.join()
