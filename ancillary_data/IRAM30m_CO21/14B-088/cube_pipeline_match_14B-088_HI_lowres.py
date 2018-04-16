
'''
Make a version of the CO(2-1) data at the resolution of 2x and 5x of the
14B-088 HI data.
'''

import numpy as np
import skimage.morphology as mo
from spectral_cube import SpectralCube, Projection
from astropy.io import fits
import os
from reproject import reproject_interp

from cube_analysis import run_pipeline
from cube_analysis.reprojection import reproject_cube

from paths import (fourteenB_HI_data_wGBT_path,
                   iram_co21_data_path,
                   iram_co21_14B088_data_path)


two_beam_hi_cubename = fourteenB_HI_data_wGBT_path("smooth_2beam/M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.fits")
vla_cube = SpectralCube.read(two_beam_hi_cubename)

# Mask out the really noisy edges so these don't get included when smoothing
co21_noise = Projection.from_hdu(fits.open(iram_co21_data_path("m33.rms.fits"))[0])
# NaNs are negative values in the map
co21_noise[co21_noise.value < 0] = np.NaN
# Initial mask below 0.04 K
mask = co21_noise.value < 0.04
# Remove holes
mask = mo.remove_small_holes(mask)
# Erode the edges by 13 pixels (beam width)
mask = mo.erosion(mask, mo.square(13))
co21_noise[~mask] = np.NaN

large_beam = vla_cube.beam

# Save to its own folder
out_folder = iram_co21_14B088_data_path("smooth_2beam", no_check=True)
if not os.path.exists(out_folder):
    os.mkdir(out_folder)

# Now smooth the noise map
co21_noise_smoothed = co21_noise.convolve_to(large_beam)
co21_noise_smoothed = co21_noise_smoothed.reproject(vla_cube[0].header)

reproj_mask = reproject_interp((mask, co21_noise.header),
                               co21_noise_smoothed.header)[0]
reproj_mask = ~(reproj_mask > 0)

# Dilate a bit to get the lowered edges from convolving
# About half of beam width
reproj_mask = mo.dilation(reproj_mask, mo.square(7))

co21_noise_smoothed[reproj_mask] = np.NaN

smooth_noise_name = os.path.join(out_folder,
                                 "m33.rms.14B-088_HI.38arcsec.fits")
co21_noise_smoothed.write(smooth_noise_name, overwrite=True)

reproject_cube(iram_co21_data_path("m33.co21_iram.fits"),
               two_beam_hi_cubename,
               "m33.co21_iram.14B-088_HI.38arcsec.fits",
               output_folder=out_folder,
               save_spectral=False,
               is_huge=False,
               reproject_type='spatial',
               common_beam=True,
               verbose=True,
               chunk=20)

run_pipeline(os.path.join(out_folder,
                          "m33.co21_iram.14B-088_HI.38arcsec.fits"),
             out_folder,
             masking_kwargs={"method": "ppv_dilation",
                             "save_cube": True,
                             "noise_map": co21_noise_smoothed.value,
                             "min_sig": 1.,
                             "max_sig": 2.5,
                             # the spatial cell of both grids are the same size
                             "min_pix": 27,
                             "verbose": True
                             },
             moment_kwargs={"num_cores": 1,
                            "verbose": True})

del vla_cube

# Now smooth to 5x the beam
five_beam_hi_cubename = fourteenB_HI_data_wGBT_path("smooth_5beam/M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.fits")
vla_cube = SpectralCube.read(five_beam_hi_cubename)

# Save to its own folder
out_folder = iram_co21_14B088_data_path("smooth_5beam", no_check=True)
if not os.path.exists(out_folder):
    os.mkdir(out_folder)

co21_noise = Projection.from_hdu(fits.open(iram_co21_data_path("m33.rms.fits"))[0])
# NaNs are negative values in the map
co21_noise[co21_noise.value < 0] = np.NaN
# Initial mask below 0.04 K
mask = co21_noise.value < 0.04

# Initial mask below 0.04 K
mask = co21_noise.value < 0.04
# Remove holes
mask = mo.remove_small_holes(mask)
# Erode the edges by 13 pixels (beam width)
mask = mo.erosion(mask, mo.square(32))
co21_noise[~mask] = np.NaN

large_beam = vla_cube.beam

# Now smooth the noise map
co21_noise_smoothed = co21_noise.convolve_to(large_beam)
co21_noise_smoothed = co21_noise_smoothed.reproject(vla_cube[0].header)

reproj_mask = reproject_interp((mask, co21_noise.header),
                               co21_noise_smoothed.header)[0]
reproj_mask = reproj_mask > 0

# Dilate a bit to get the lowered edges from convolving
# About half of beam width
reproj_mask = mo.dilation(reproj_mask, mo.square(16))

co21_noise_smoothed[~reproj_mask] = np.NaN

smooth_noise_name = os.path.join(out_folder,
                                 "m33.rms.14B-088_HI.95arcsec.fits")
co21_noise_smoothed.write(smooth_noise_name, overwrite=True)

reproject_cube(iram_co21_data_path("m33.co21_iram.fits"),
               five_beam_hi_cubename,
               "m33.co21_iram.14B-088_HI.95arcsec.fits",
               output_folder=out_folder,
               save_spectral=False,
               is_huge=True,
               reproject_type='spatial',
               common_beam=True,
               verbose=True,
               chunk=20)

# Use the noise mask to mask the reprojected cube
cube = SpectralCube.read(os.path.join(out_folder, "m33.co21_iram.14B-088_HI.95arcsec.fits"))
cube = cube.with_mask(reproj_mask)
cube.write(os.path.join(out_folder, "m33.co21_iram.14B-088_HI.95arcsec.fits"),
           overwrite=True)
del cube

run_pipeline(os.path.join(out_folder,
                          "m33.co21_iram.14B-088_HI.95arcsec.fits"),
             out_folder,
             masking_kwargs={"method": "ppv_dilation",
                             "save_cube": True,
                             "noise_map": co21_noise_smoothed.value,
                             "min_sig": 0.5,
                             "max_sig": 2.,
                             # the spatial cell of both grids are the same size
                             "min_pix": 27,
                             "verbose": True
                             },
             moment_kwargs={"num_cores": 1,
                            "verbose": True})
