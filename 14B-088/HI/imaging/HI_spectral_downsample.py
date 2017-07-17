
'''
Create a down-sampled version of the HI cube.
'''

import os
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import largest_beam
from astropy import log
from astropy.convolution import Box1DKernel
import numpy as np

from cube_analysis.io_utils import save_to_huge_fits

from paths import (fourteenB_HI_data_wGBT_path,
                   fourteenB_wGBT_HI_file_dict)


cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])

# Smooth to the largest beam
log.info("Convolving to largest beam.")
cube = cube.convolve_to(largest_beam(cube.beams))

# Smooth over 5 channels to get an effective channel width of ~ 1 km/s
num_chans = 5
chan_width = np.diff(cube.spectral_axis[:2])[0]
kernel = Box1DKernel(num_chans)

# Smooth first
log.info("Spectrally smoothing.")
new_cube = cube.spectral_smooth(kernel)

# Then interpolate
log.info("Spectral interpolating.")
new_spec_length = cube.shape[0] / num_chans if cube.shape[0] % num_chans == 0 \
    else (cube.shape[0] / num_chans) + 1
new_specaxis = np.linspace(cube.spectral_axis[0].value,
                           cube.spectral_axis[-1].value,
                           new_spec_length,
                           endpoint=True) * new_cube.spectral_axis.unit
new_cube = new_cube.spectral_interpolate(new_specaxis)

log.info("Writing out.")
out_folder = fourteenB_HI_data_wGBT_path("downsamp_1kms", no_check=True)
if not os.path.exists(out_folder):
    os.mkdir(out_folder)
save_to_huge_fits(fourteenB_HI_data_wGBT_path('downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.1kms.fits',
                                              no_check=True),
                  new_cube, chunk=50)

# Do the same for the rotation-corrected cube
cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['RotSub_Cube'])

# Smooth over 5 channels to get an effective channel width of ~ 1 km/s
num_chans = 5
chan_width = np.diff(cube.spectral_axis[:2])[0]
kernel = Box1DKernel(num_chans)

# Smooth first
log.info("Spectrally smoothing.")
new_cube = cube.spectral_smooth(kernel)

# Then interpolate
log.info("Spectral interpolating.")
new_spec_length = cube.shape[0] / num_chans if cube.shape[0] % num_chans == 0 \
    else (cube.shape[0] / num_chans) + 1
new_specaxis = np.linspace(cube.spectral_axis[0].value,
                           cube.spectral_axis[-1].value,
                           new_spec_length,
                           endpoint=True) * new_cube.spectral_axis.unit
new_cube = new_cube.spectral_interpolate(new_specaxis)

log.info("Writing out.")
out_folder = fourteenB_HI_data_wGBT_path("downsamp_1kms", no_check=True)
if not os.path.exists(out_folder):
    os.mkdir(out_folder)
save_to_huge_fits(fourteenB_HI_data_wGBT_path('downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.rotation_corrected.1kms.fits',
                                              no_check=True),
                  new_cube, chunk=50)
