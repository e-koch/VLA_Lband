
'''
Convolve the VLA + GBT data to 2 * beam and 5 * beam, then run the
masking and moments pipeline.

Make signal masks and compute the moments.
'''

from astropy import log
import os
from radio_beam import Beam
from spectral_cube import SpectralCube
from cube_analysis import run_pipeline

from paths import (fourteenB_wGBT_HI_file_dict, fourteenB_HI_data_wGBT_path)

file_path = fourteenB_HI_data_wGBT_path("smooth_2beam", no_check=True)
if not os.path.exists(file_path):
    os.mkdir(file_path)

cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["Cube"])
# Convolve to 2 * beam. May as well make it circular.
beam2 = Beam(2 * cube.beams.largest_beam().major)

conv_cube = cube.convolve_to(beam2)

file_name = os.path.join(file_path, "M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.fits")
conv_cube.write(file_name)

del conv_cube
del cube

log.info("Masking and moments for the VLA+GBT 2 * beam cube")
run_pipeline(file_name,
             file_path,
             masking_kwargs={"method": "ppv_connectivity",
                             "save_cube": True,
                             "is_huge": True,
                             "noise_map": None,
                             "smooth_chans": 31,
                             "min_chan": 10,
                             "peak_snr": 5.,
                             "min_snr": 2,
                             "edge_thresh": 1,
                             "verbose": False,
                             "show_plots": False,
                             },
             moment_kwargs={"num_cores": 6,
                            "verbose": True})

print("Running 5 beam convolution and masking.")

file_path = fourteenB_HI_data_wGBT_path("smooth_5beam", no_check=True)
if not os.path.exists(file_path):
    os.mkdir(file_path)

cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["Cube"])

beam5 = Beam(5 * cube.beams.largest_beam().major)

conv_cube = cube.convolve_to(beam5)

file_name = os.path.join(file_path, "M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.fits")
conv_cube.write(file_name)

del conv_cube
del cube

# VLA+GBT cube
log.info("Masking and moments for the VLA+GBT 5 * beam cube")
run_pipeline(file_name,
             file_path,
             masking_kwargs={"method": "ppv_connectivity",
                             "save_cube": True,
                             "is_huge": True,
                             "noise_map": None,
                             "smooth_chans": 31,
                             "min_chan": 10,
                             "peak_snr": 5.,
                             "min_snr": 2,
                             "edge_thresh": 1,
                             },
             moment_kwargs={"num_cores": 1,
                            "verbose": True,
                            "chunk_size": 2e5})
