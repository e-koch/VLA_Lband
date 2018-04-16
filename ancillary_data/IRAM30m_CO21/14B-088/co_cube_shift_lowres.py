
'''
Shift the CO data with different velocity surfaces.
'''

from spectral_cube import SpectralCube
from astropy.io import fits
from astropy import log
import astropy.units as u
import numpy as np
import scipy.ndimage as nd
from os.path import join as osjoin

from cube_analysis.spectra_shifter import cube_shifter

from paths import (iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path)

from galaxy_params import gal_feath as gal


num_cores = 4
chunk_size = 50000
verbose = False

smooth_2beam_path_co = lambda x: osjoin(iram_co21_14B088_data_path("smooth_2beam"), x)
smooth_2beam_path_hi = lambda x: osjoin(fourteenB_HI_data_wGBT_path("smooth_2beam"), x)

cube = SpectralCube.read(smooth_2beam_path_co("m33.co21_iram.14B-088_HI.38arcsec.fits"))
mask = SpectralCube.read(smooth_2beam_path_co("m33.co21_iram.14B-088_HI.38arcsec_source_mask.fits"))
rms_map = fits.open(smooth_2beam_path_co("m33.rms.14B-088_HI.38arcsec.fits"))[0]

rms_mask = np.isfinite(rms_map.data)

# Reduce the edges to make sure we are avoiding the severely noisy regions
valid_mask = nd.binary_erosion(rms_mask, np.ones((3, 3)), iterations=10)

mom1 = fits.open(smooth_2beam_path_hi("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.mom1.fits"))[0].data * u.m / u.s
mom1[~valid_mask] = np.NaN
peakvels = \
    fits.open(smooth_2beam_path_hi("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.peakvels.fits"))[0].data * u.m / u.s
peakvels[~valid_mask] = np.NaN

# Shift with the first moment
log.info("Shifting with HI centroid.")
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=smooth_2beam_path_co("m33.co21_iram.14B-088_HI.38arcsec.centroid_corrected.fits"),
             return_spectra=False, verbose=verbose,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI centroid.")
cube_shifter(mask, mom1, gal.vsys, save_shifted=True,
             save_name=smooth_2beam_path_co("m33.co21_iram.14B-088_HI.38arcsec_source_mask.centroid_corrected.fits"),
             return_spectra=False, verbose=verbose, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Shift with the HI peak velocities
log.info("Shifting with HI peak velocity.")
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=smooth_2beam_path_co("m33.co21_iram.14B-088_HI.38arcsec.peakvels_corrected.fits"),
             return_spectra=False, verbose=verbose,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI peak velocity.")
cube_shifter(mask, peakvels, gal.vsys, save_shifted=True,
             save_name=smooth_2beam_path_co("m33.co21_iram.14B-088_HI.38arcsec_source_mask.peakvels_corrected.fits"),
             return_spectra=False, verbose=verbose, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

del cube, mask, rms_map

# Now at 5 HI beam

smooth_5beam_path_co = lambda x: osjoin(iram_co21_14B088_data_path("smooth_5beam"), x)
smooth_5beam_path_hi = lambda x: osjoin(fourteenB_HI_data_wGBT_path("smooth_5beam"), x)

cube = SpectralCube.read(smooth_5beam_path_co("m33.co21_iram.14B-088_HI.95arcsec.fits"))
mask = SpectralCube.read(smooth_5beam_path_co("m33.co21_iram.14B-088_HI.95arcsec_source_mask.fits"))
rms_map = fits.open(smooth_5beam_path_co("m33.rms.14B-088_HI.95arcsec.fits"))[0]

rms_mask = np.isfinite(rms_map.data)

# Reduce the edges to make sure we are avoiding the severely noisy regions
valid_mask = nd.binary_erosion(rms_mask, np.ones((3, 3)), iterations=10)

# Make the same versions using the feathered HI data
mom1 = fits.open(smooth_5beam_path_hi("M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.mom1.fits"))[0].data * u.m / u.s
mom1[~valid_mask] = np.NaN
peakvels = \
    fits.open(smooth_5beam_path_hi("M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.peakvels.fits"))[0].data * u.m / u.s
peakvels[~valid_mask] = np.NaN

# Shift with the first moment
log.info("Shifting with HI centroid. 5 * beam")
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=smooth_5beam_path_co("m33.co21_iram.14B-088_HI.95arcsec.centroid_corrected.fits"),
             return_spectra=False, verbose=verbose,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI centroid. 5 * beam")
cube_shifter(mask, mom1, gal.vsys, save_shifted=True,
             save_name=smooth_5beam_path_co("m33.co21_iram.14B-088_HI.95arcsec_source_mask.centroid_corrected.fits"),
             return_spectra=False, verbose=verbose, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Shift with the HI peak velocities
log.info("Shifting with HI peak velocity. 5 * beam")
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=smooth_5beam_path_co("m33.co21_iram.14B-088_HI.95arcsec.peakvels_corrected.fits"),
             return_spectra=False, verbose=verbose,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI peak velocity. 5 * beam")
cube_shifter(mask, peakvels, gal.vsys, save_shifted=True,
             save_name=smooth_5beam_path_co("m33.co21_iram.14B-088_HI.95arcsec_source_mask.peakvels_corrected.fits"),
             return_spectra=False, verbose=verbose, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)
