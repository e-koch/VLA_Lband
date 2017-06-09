
'''
Create feathered VLA cubes with the 3 different SD datasets.
'''

from spectral_cube import SpectralCube
import os
from astropy import log

from cube_analysis.feather_cubes import feather_cube

from paths import fourteenB_HI_data_path, data_path
from constants import hi_freq

# Load the non-pb masked cube
vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))

# Set which of the cubes to feather
run_gbt_highres = True
run_gbt = True
run_ebhis = True
run_arecibo = True

num_cores = 6

if run_gbt_highres:
    log.info("Feathering with high res. GBT")

    gbt_path = os.path.join(data_path, "GBT")
    gbt_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid_registered.fits")
    gbt_cube = SpectralCube.read(gbt_name)

    output_path = os.path.join(data_path, "VLA/14B-088/HI/full_imaging_wGBT/")
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    save_name = os.path.join(output_path, "M33_14B-088_HI.clean.image.GBT_feathered.fits")

    feather_cube(vla_cube, gbt_cube, frequency=hi_freq, save_feather=True,
                 save_name=save_name, num_cores=num_cores)


if run_gbt:
    log.info("Feathering with low res. GBT")

    gbt_path = os.path.join(data_path, "GBT")
    gbt_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_Tmb_14B088_spectralregrid_registered.fits")
    gbt_cube = SpectralCube.read(gbt_name)

    output_path = os.path.join(data_path, "VLA/14B-088/HI/full_imaging_wGBT_lowres/")
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    save_name = os.path.join(output_path, "M33_14B-088_HI.clean.image.GBT_lowres_feathered.fits")

    feather_cube(vla_cube, gbt_cube, frequency=hi_freq, save_feather=True,
                 save_name=save_name, num_cores=num_cores)

if run_ebhis:
    log.info("Feathering with EBHIS")

    ebhis_path = os.path.join(data_path, "EBHIS")
    ebhis_name = os.path.join(ebhis_path, "14B-088_items/m33_ebhis_14B088_spectralregrid_registered.fits")
    ebhis_cube = SpectralCube.read(ebhis_name)

    output_path = os.path.join(data_path, "VLA/14B-088/HI/full_imaging_wEBHIS/")
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    save_name = os.path.join(output_path, "M33_14B-088_HI.clean.image.EBHIS_feathered.fits")

    feather_cube(vla_cube, ebhis_cube, frequency=hi_freq, save_feather=True,
                 save_name=save_name, num_cores=num_cores)

if run_arecibo:
    log.info("Feathering with Arecibo")

    arecibo_path = os.path.join(data_path, "Arecibo")
    arecibo_name = os.path.join(arecibo_path, "14B-088_items_new/m33_arecibo_14B088_spectralregrid_registered.fits")
    arecibo_cube = SpectralCube.read(arecibo_name)

    output_path = os.path.join(data_path, "VLA/14B-088/HI/full_imaging_wArecibo/")
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    save_name = os.path.join(output_path, "M33_14B-088_HI.clean.image.Arecibo_feathered.fits")

    feather_cube(vla_cube, arecibo_cube, frequency=hi_freq, save_feather=True,
                 save_name=save_name, num_cores=num_cores)
