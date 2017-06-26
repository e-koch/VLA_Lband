
'''
Make a version of the CO(2-1) data at the resolution of the 14B-088 HI data
'''

import numpy as np
import skimage.morphology as mo
from spectral_cube import SpectralCube, Projection
from spectral_cube.cube_utils import largest_beam
from astropy.io import fits
import os

from cube_analysis import run_pipeline
from cube_analysis.reprojection import reproject_cube

from paths import fourteenB_HI_file_dict, iram_co21_data_path


vla_cube = SpectralCube.read(fourteenB_HI_file_dict["Cube"])

# co21_cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))

# Mask out the really noisy edges so these don't get included when smoothing
co21_noise = Projection.from_hdu(fits.open(iram_co21_data_path("m33.rms.fits"))[0])
# NaNs are negative values in the map
co21_noise[co21_noise.value < 0] = np.NaN
# Initial mask below 0.04 K
mask = co21_noise.value < 0.04
# Remove holes
mask = mo.remove_small_holes(mask)
# Erode the edges by 10 pixels (chosen by-eye)
mask = mo.erosion(mask, mo.square(10))
co21_noise[~mask] = np.NaN

# Save this version of the noise
co21_noise.write(iram_co21_data_path("m33.rms.masked.fits", no_check=True),
                 overwrite=True)

# Smooth to the largest beam (they differ by tiny fractions anyway)
large_beam = largest_beam(vla_cube.beams)

# Save to its own folder
out_folder = iram_co21_data_path("14B-088", no_check=True)
if not os.path.exists(out_folder):
    os.mkdir(out_folder)

# Now smooth the noise map
co21_noise_smoothed = co21_noise.convolve_to(large_beam)
co21_noise_smoothed = co21_noise_smoothed.reproject(vla_cube[0].header)
reproj_mask = co21_noise_smoothed < np.nanmin(co21_noise) / 2.
# Dilate a bit to get the lowered edges from convolving
reproj_mask = mo.dilation(reproj_mask, mo.square(7))

co21_noise_smoothed[reproj_mask] = np.NaN

smooth_noise_name = iram_co21_data_path("14B-088/m33.rms.14B-088_HI.fits",
                                        no_check=True)
co21_noise_smoothed.write(smooth_noise_name, overwrite=True)

reproject_cube(iram_co21_data_path("m33.co21_iram.fits"),
               fourteenB_HI_file_dict["Cube"],
               "m33.co21_iram.14B-088_HI.fits",
               output_folder=out_folder,
               save_spectral=False,
               is_huge=True,
               reproject_type='spatial',
               common_beam=True,
               verbose=True,
               chunk=20)

run_pipeline(os.path.join(out_folder, "m33.co21_iram.14B-088_HI.fits"),
             out_folder,
             masking_kwargs={"method": "ppv_dilation",
                             "save_cube": True,
                             "noise_map": co21_noise_smoothed.value,
                             "min_sig": 2.2,
                             "max_sig": 4,
                             # the spatial cell of both grids are the same size
                             "min_pix": 27,
                             "verbose": True
                             },
             moment_kwargs={"num_cores": 6,
                            "verbose": True})

# We'll also make a smoothed, non-reprojected version.
co_cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
co_cube = co_cube.convolve_to(large_beam)
co_cube.write(iram_co21_data_path("14B-088/m33.co21_iram.14B-088_HI.smoothonly.fits"),
              overwrite=True)
