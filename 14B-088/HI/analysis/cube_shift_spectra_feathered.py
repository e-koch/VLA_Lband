
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
from multiprocessing import Pool

from cube_analysis.spectra_shifter import cube_shifter

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict)
from galaxy_params import gal_feath as gal


cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])

mom1 = fits.open(fourteenB_wGBT_HI_file_dict["Moment1"])[0].data * u.m / u.s
peakvels = \
    fits.open(fourteenB_wGBT_HI_file_dict["PeakVels"])[0].data * u.m / u.s
rotmod_name = fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_noradial_nowarp_output/rad.fitmod.fits")
rot_model = fits.open(rotmod_name)[0].data * u.m / u.s

pool = Pool(6, maxtasksperchild=5000)

# Shift with the first moment
log.info("Shifting w/ centroid")
centroidsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.centroid_corrected.fits"
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(centroidsub_cube_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, pool=pool)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ centroid")
mask_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["Source_Mask"])
centroidsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked_source_mask.centroid_corrected.fits"
cube_shifter(mask_cube, mom1, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(centroidsub_mask_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, pool=pool,
             is_mask=True)

# Shift wrt. peak velocities
log.info("Shifting w/ peak velocity")
peakvelsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peakvels_corrected.fits"
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(peakvelsub_cube_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, pool=pool)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ peak velocity")
peakvelsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked_source_mask.peakvels_corrected.fits"
cube_shifter(mask_cube, peakvels, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(peakvelsub_mask_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, pool=pool,
             is_mask=True)

# Shift with the rotation model
log.info("Shifting w/ rotation velocity")
rotsub_cube_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.rotation_corrected.fits"
cube_shifter(cube, rot_model, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(rotsub_cube_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, pool=pool)

# Also shift the signal mask to match those shifted here.
log.info("Shifting mask w/ rotation velocity")
rotsub_mask_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked_source_mask.rotation_corrected.fits"
cube_shifter(mask_cube, rot_model, gal.vsys, save_shifted=True,
             save_name=fourteenB_HI_data_wGBT_path(rotsub_mask_name,
                                                   no_check=True),
             return_spectra=False, verbose=True, pool=pool,
             is_mask=True)
pool.close()
pool.join()
