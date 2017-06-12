
'''
Make a version of the CO(2-1) data at the resolution of the 14B-088 HI data
'''

from spectral_cube import SpectralCube, Projection
from spectral_cube.cube_utils import largest_beam
from astropy.io import fits
import os

from cube_analysis import run_pipeline

from paths import fourteenB_HI_file_dict, iram_co21_data_path


vla_cube = SpectralCube.read(fourteenB_HI_file_dict["Cube"])

co21_cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))

# Smooth to the largest beam (they differ by tiny fractions anyway)

large_beam = largest_beam(vla_cube.beams)
# smoothed_cube = co21_cube.convolve_to(large_beam)

# Save to its own folder
if not os.path.exists(iram_co21_data_path("14B-088", no_check=True)):
    os.mkdir(iram_co21_data_path("14B-088", no_check=True))

smooth_name = iram_co21_data_path("14B-088/m33.co21_iram.14B-088_HI_smoothed.fits", no_check=True)
# smoothed_cube.write(smooth_name)

# Now smooth the noise map
co21_noise = Projection.from_hdu(fits.open(iram_co21_data_path("m33.rms.fits"))[0])
# co21_noise_smoothed = co21_noise.convolve_to(large_beam)
smooth_noise_name = iram_co21_data_path("14B-088/m33.rms.14B-088_HI_smoothed.fits", no_check=True)
# co21_noise_smoothed.write(smooth_noise_name)

# Find a signal mask and derive moment arrays
run_pipeline(smooth_name, iram_co21_data_path("14B-088"),
             masking_kwargs={"method": "ppv_dilation",
                             "save_cube": True,
                             "noise_map": smooth_noise_name,
                             "min_sig": 3,
                             "max_sig": 5,
                             "min_pix": 27,
                             },
             moment_kwargs={"num_cores": 6,
                            "verbose": True})
