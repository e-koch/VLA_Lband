
'''
Shift the low-resolution versions of the combined HI cube.
Shifting w.r.t. the peak and centroid

'''

from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
from astropy import log
from os.path import join as osjoin

from cube_analysis.spectra_shifter import cube_shifter

from paths import (fourteenB_HI_data_wGBT_path)
from galaxy_params import gal_feath as gal


num_cores = 4

# Set the number of spectra to load in and operate on at once.
chunk_size = 200000

# 2 * beam = 38'' cube
smooth_2beam_folder = fourteenB_HI_data_wGBT_path("smooth_2beam")

cube = SpectralCube.read(osjoin(smooth_2beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.fits"))
mask_cube = SpectralCube.read(osjoin(smooth_2beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec_source_mask.fits"))

mom1 = fits.open(osjoin(smooth_2beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.mom1.fits"))[0].data * u.m / u.s
peakvels = \
    fits.open(osjoin(smooth_2beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.peakvels.fits"))[0].data * u.m / u.s

# Shift with the first moment
log.info("Shifting w/ centroid")
centroidsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.centroid_corrected.fits"
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_2beam_folder, centroidsub_cube_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             chunk_size=chunk_size, pad_edges=True)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ centroid")
centroidsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec_source_mask.centroid_corrected.fits"
cube_shifter(mask_cube, mom1, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_2beam_folder, centroidsub_mask_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             is_mask=True, chunk_size=chunk_size, pad_edges=True)

# Shift wrt. peak velocities
log.info("Shifting w/ peak velocity")
peakvelsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.peakvels_corrected.fits"
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_2beam_folder, peakvelsub_cube_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             chunk_size=chunk_size, pad_edges=True)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ peak velocity")
peakvelsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec_source_mask.peakvels_corrected.fits"
cube_shifter(mask_cube, peakvels, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_2beam_folder, peakvelsub_mask_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             is_mask=True, chunk_size=chunk_size, pad_edges=True)

del cube
del mask_cube

# 5 * beam = 95'' cube
smooth_5beam_folder = fourteenB_HI_data_wGBT_path("smooth_5beam")

cube = SpectralCube.read(osjoin(smooth_5beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.fits"))
mask_cube = SpectralCube.read(osjoin(smooth_5beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec_source_mask.fits"))

mom1 = fits.open(osjoin(smooth_5beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.mom1.fits"))[0].data * u.m / u.s
peakvels = \
    fits.open(osjoin(smooth_5beam_folder, "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.peakvels.fits"))[0].data * u.m / u.s

# Shift with the first moment
log.info("Shifting w/ centroid")
centroidsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.centroid_corrected.fits"
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_5beam_folder, centroidsub_cube_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             chunk_size=chunk_size, pad_edges=True)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ centroid")
centroidsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec_source_mask.centroid_corrected.fits"
cube_shifter(mask_cube, mom1, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_5beam_folder, centroidsub_mask_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             is_mask=True, chunk_size=chunk_size, pad_edges=True)

# Shift wrt. peak velocities
log.info("Shifting w/ peak velocity")
peakvelsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.peakvels_corrected.fits"
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_5beam_folder, peakvelsub_cube_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             chunk_size=chunk_size, pad_edges=True)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ peak velocity")
peakvelsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec_source_mask.peakvels_corrected.fits"
cube_shifter(mask_cube, peakvels, gal.vsys, save_shifted=True,
             save_name=osjoin(smooth_5beam_folder, peakvelsub_mask_name),
             return_spectra=False, verbose=True, num_cores=num_cores,
             is_mask=True, chunk_size=chunk_size, pad_edges=True)
