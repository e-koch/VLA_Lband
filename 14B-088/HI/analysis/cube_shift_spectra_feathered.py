
'''
Given a velocity surface, shift all spectra in the feathered cube.

Mask versions shifted w.r.t:

 * centroids
 * peak velocities
 * rotation model

'''

from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
from astropy import log
import numpy as np

from cube_analysis.spectra_shifter import cube_shifter

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict)
from galaxy_params import gal_feath as gal


num_cores = 4

# Set the number of spectra to load in and operate on at once.
chunk_size = 50000

cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])
mask_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["Source_Mask"])

mom1 = fits.open(fourteenB_wGBT_HI_file_dict["Moment1"])[0].data * u.m / u.s
peakvels = \
    fits.open(fourteenB_wGBT_HI_file_dict["PeakVels"])[0].data * u.m / u.s
rotmod_name = fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_noradial_nowarp_output/rad.fitmod.fits")
rot_model = fits.open(rotmod_name)[0].data * u.m / u.s
# Don't bother using points outside of the cube mask.
rot_model[np.isnan(peakvels)] = np.NaN

# Shift with the first moment
log.info("Shifting w/ centroid")
centroidsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.centroid_corrected.fits"
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(centroidsub_cube_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, num_cores=num_cores,
             chunk_size=chunk_size)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ centroid")
centroidsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked_source_mask.centroid_corrected.fits"
cube_shifter(mask_cube, mom1, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(centroidsub_mask_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, num_cores=num_cores,
             is_mask=True, chunk_size=chunk_size)

# Shift wrt. peak velocities
log.info("Shifting w/ peak velocity")
peakvelsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peakvels_corrected.fits"
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(peakvelsub_cube_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, num_cores=num_cores,
             chunk_size=chunk_size)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ peak velocity")
peakvelsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked_source_mask.peakvels_corrected.fits"
cube_shifter(mask_cube, peakvels, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(peakvelsub_mask_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, num_cores=num_cores,
             is_mask=True, chunk_size=chunk_size)

# Shift with the rotation model
log.info("Shifting w/ rotation velocity")
rotsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.rotation_corrected.fits"
cube_shifter(cube, rot_model, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(rotsub_cube_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, num_cores=num_cores,
             chunk_size=chunk_size)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ rotation velocity")
rotsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked_source_mask.rotation_corrected.fits"
cube_shifter(mask_cube, rot_model, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(rotsub_mask_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, num_cores=num_cores,
             is_mask=True, chunk_size=chunk_size)
