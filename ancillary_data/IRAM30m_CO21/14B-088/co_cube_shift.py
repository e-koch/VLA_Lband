
'''
Shift the CO data with different velocity surfaces.
'''

from spectral_cube import SpectralCube
from astropy.io import fits
from astropy import log
import astropy.units as u
import numpy as np
import sys

from cube_analysis.spectra_shifter import cube_shifter

from paths import (iram_co21_data_path, iram_co21_14B088_reproj_data_path,
                   fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   fourteenB_HI_data_path, fourteenB_HI_data_wGBT_path,
                   c_hi_analysispath)

sys.path.append(c_hi_analysispath(""))
from galaxy_params import gal_feath as gal


num_cores = 4
chunk_size = 50000

cube = SpectralCube.read(iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.fits"))
mask = SpectralCube.read(iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_source_mask.fits"))

valid_mask = mask.sum(0).value > 0

mom1 = fits.open(fourteenB_HI_file_dict["Moment1"])[0].data * u.m / u.s
mom1[~valid_mask] = np.NaN
peakvels = \
    fits.open(fourteenB_HI_file_dict["PeakVels"])[0].data * u.m / u.s
peakvels[~valid_mask] = np.NaN
rotmod_name = fourteenB_HI_data_path("diskfit_peakvels_noasymm_noradial_nowarp_output/rad.fitmod.fits")
rot_model = fits.open(rotmod_name)[0].data * u.m / u.s
# Don't bother using points outside of the cube mask.
rot_model[~valid_mask] = np.NaN

log.info("Shifting with the VLA-only velocity maps.")

# Shift with the first moment
log.info("Shifting with HI centroid.")
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.centroid_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI centroid.")
cube_shifter(mask, mom1, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_source_mask.centroid_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Shift with the HI peak velocities
log.info("Shifting with HI peak velocity.")
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.peakvels_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI peak velocity.")
cube_shifter(mask, peakvels, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_source_mask.peakvels_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Shift with the rotation curve
log.info("Shifting with HI rotation model.")
cube_shifter(cube, rot_model, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.rotation_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI rotation model.")
cube_shifter(mask, rot_model, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_source_mask.rotation_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Make the same versions using the feathered HI data
mom1 = fits.open(fourteenB_wGBT_HI_file_dict["Moment1"])[0].data * u.m / u.s
mom1[~valid_mask] = np.NaN
peakvels = \
    fits.open(fourteenB_wGBT_HI_file_dict["PeakVels"])[0].data * u.m / u.s
peakvels[~valid_mask] = np.NaN
rotmod_name = fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_noradial_nowarp_output/rad.fitmod.fits")
rot_model = fits.open(rotmod_name)[0].data * u.m / u.s
# Don't bother using points outside of the cube mask.
rot_model[~valid_mask] = np.NaN

log.info("Shifting with the feathered velocity maps.")

# Shift with the first moment
log.info("Shifting with HI centroid.")
cube_shifter(cube, mom1, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_feather.centroid_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI centroid.")
cube_shifter(mask, mom1, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_feather_source_mask.centroid_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Shift with the HI peak velocities
log.info("Shifting with HI peak velocity.")
cube_shifter(cube, peakvels, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_feather.peakvels_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI peak velocity.")
cube_shifter(mask, peakvels, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_feather_source_mask.peakvels_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Shift with the rotation curve
log.info("Shifting with HI rotation model.")
cube_shifter(cube, rot_model, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_feather.rotation_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)
log.info("Shifting mask with HI rotation model.")
cube_shifter(mask, rot_model, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_feather_source_mask.rotation_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Shift with the CO peak velocities & centroid
# Do this for the different CO cubes

# Smoothed and reprojected version
log.info("Shifting the reprojected cube now.")
peakvels_co = fits.open(iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.peakvels.fits"))[0].data * u.m / u.s
centroid_co = fits.open(iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.mom1.fits"))[0].data * u.m / u.s

log.info("Shifting with the CO peak velocity.")
cube_shifter(cube, peakvels_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.co_peakvels_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)

log.info("Shifting mask with the CO peak velocity.")
cube_shifter(mask, peakvels_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_source_mask.co_peakvels_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

log.info("Shifting with the CO centroid.")
cube_shifter(cube, centroid_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj.co_centroid_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)

log.info("Shifting mask with the CO centroid.")
cube_shifter(mask, centroid_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_14B088_reproj_data_path("m33.co21_iram.14B-088_HI_reproj_source_mask.centroid_corrected.fits",
                                                         no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

# Run on the original CO cube
cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
mask = SpectralCube.read(iram_co21_data_path("m33.co21_iram_source_mask.fits"))

peakvels_co = fits.open(iram_co21_data_path("m33.co21_iram.peakvels.fits"))[0].data * u.m / u.s
centroid_co = fits.open(iram_co21_data_path("m33.co21_iram.mom1.fits"))[0].data * u.m / u.s

log.info("Shifting the original cube version.")
log.info("Shifting with the CO peak velocity.")
cube_shifter(cube, peakvels_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_data_path("m33.co21_iram.co_peakvels_corrected.fits",
                                           no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)

log.info("Shifting mask with the CO peak velocity.")
cube_shifter(mask, peakvels_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_data_path("m33.co21_iram_source_mask.co_peakvels_corrected.fits",
                                           no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)

log.info("Shifting with the CO centroid.")
cube_shifter(cube, centroid_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_data_path("m33.co21_iram.co_centroid_corrected.fits",
                                           no_check=True),
             return_spectra=False, verbose=True,
             num_cores=num_cores, chunk_size=chunk_size)

log.info("Shifting mask with the CO centroid.")
cube_shifter(mask, centroid_co, gal.vsys, save_shifted=True,
             save_name=iram_co21_data_path("m33.co21_iram_source_mask.co_centroid_corrected.fits",
                                           no_check=True),
             return_spectra=False, verbose=True, is_mask=True,
             num_cores=num_cores, chunk_size=chunk_size)
