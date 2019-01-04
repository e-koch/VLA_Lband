
'''
Match the GBT cube to the 14B footprint and feather with the 5 km/s channel
cube.
'''

import os
from os.path import join as osjoin
from spectral_cube import SpectralCube

from cube_analysis.feather_cubes import feather_cube

from paths import data_path
from constants import hi_freq


# 5 km/s path
fivekms_noSD_path = osjoin(data_path, "VLA/14B-088/HI/full_imaging_5kms_noSD")

vla_cube = SpectralCube.read(osjoin(fivekms_noSD_path, "M33_14B-088_HI_5kms.clean.image.pbcor.fits"))

num_cores = 4

gbt_path = os.path.join(data_path, "GBT")
gbt_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid_registered.fits")
gbt_cube = SpectralCube.read(gbt_name)

# First regrid the velocity to match the 5 km/s VLA cube
# The low spatial resolution gives smooth spectra, so don't worry about
# smoothing prior to resampling

gbt_smooth_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_5kms_registered.fits")

if not os.path.exists(gbt_smooth_name):
    gbt_smooth_cube = gbt_cube.spectral_interpolate(vla_cube.spectral_axis)
    gbt_smooth_cube.write(gbt_smooth_name)
else:
    gbt_smooth_cube = SpectralCube.read(gbt_smooth_name)

output_path = os.path.join(data_path, "VLA/14B-088/HI/full_imaging_5kms_wGBT/")
if not os.path.exists(output_path):
    os.mkdir(output_path)

save_name = os.path.join(output_path, "M33_14B-088_HI_5kms.clean.image.pbcor.GBT_feathered.fits")

feather_cube(vla_cube, gbt_smooth_cube, restfreq=hi_freq, save_feather=True,
             save_name=save_name, num_cores=num_cores)
