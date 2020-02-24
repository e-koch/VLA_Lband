
'''
Reproject onto the ACA CO(2-1) mosaic.

Since we have different versions of the full mosaic, we're
only going to make two version of the reprojected HI maps:

1) One to the ACA map without the highly asymmetric beam mosaics.
The HI map will not be convolved to a matching beam since they are already
quite similar. (But this can be checked later)

2) One to the full ACA map and convolved to its round beam of ~11 arcsec.

Also, these are just spatial reprojections, not spectral. So the CO
channel width won't matter.

'''


import os
from os.path import join as osjoin

from cube_analysis.reprojection import reproject_cube
from cube_analysis.run_pipe import run_pipeline

from paths import (seventeenB_HI_data_1kms_wGBT_path,
                   seventeenB_1kms_wGBT_HI_file_dict, aca_co21_data_path)


out_folder = seventeenB_HI_data_1kms_wGBT_path("aca_co21_match",
                                               no_check=True)

if not os.path.exists(out_folder):
    os.mkdir(out_folder)

run_noasymm = True

run_fullmos_round = True


if run_noasymm:

    out_name = seventeenB_1kms_wGBT_HI_file_dict['Cube'].split("/")[-1].rstrip(".fits") + \
        ".aca_excasymm_spatialmatch.fits"

    targ_cube = aca_co21_data_path("full_mosaic/12CO21/M33_ACA_12CO21_2p6kms_excludeasymmbeams_commonbeam.image_K.fits")

    reproject_cube(seventeenB_1kms_wGBT_HI_file_dict['Cube'],
                   targ_cube,
                   out_name,
                   output_folder=out_folder,
                   save_spectral=False,
                   is_huge=True,
                   reproject_type='spatial',
                   common_beam=False,
                   verbose=True,
                   chunk=40)

    run_pipeline(osjoin(out_folder, out_name),
                 out_folder,
                 masking_kwargs={"method": "ppv_connectivity",
                                 "save_cube": True,
                                 "is_huge": True,
                                 "smooth_chans": 6,
                                 "min_chan": 4,
                                 "peak_snr": 4.,
                                 "min_snr": 2,
                                 "edge_thresh": 1,
                                 "verbose": True,
                                 },
                 moment_kwargs={"num_cores": 1,
                                "verbose": False,
                                "chunk_size": 1e5,
                                "make_peakvels": False},)

if run_fullmos_round:

    out_name = seventeenB_1kms_wGBT_HI_file_dict['Cube'].split("/")[-1].rstrip(".fits") + \
        ".aca_fullmosaic_spatialmatch.fits"

    targ_cube = aca_co21_data_path("full_mosaic/12CO21/M33_ACA_12CO21_2p6kms_fullmosaic_roundbeam.image_K.fits")

    reproject_cube(seventeenB_1kms_wGBT_HI_file_dict['Cube'],
                   targ_cube,
                   out_name,
                   output_folder=out_folder,
                   save_spectral=False,
                   is_huge=True,
                   reproject_type='spatial',
                   common_beam=True,
                   verbose=True,
                   chunk=40)

    run_pipeline(osjoin(out_folder, out_name),
                 out_folder,
                 masking_kwargs={"method": "ppv_connectivity",
                                 "save_cube": True,
                                 "is_huge": True,
                                 "smooth_chans": 6,
                                 "min_chan": 4,
                                 "peak_snr": 4.,
                                 "min_snr": 2,
                                 "edge_thresh": 1,
                                 "verbose": True,
                                 },
                 moment_kwargs={"num_cores": 1,
                                "verbose": False,
                                "chunk_size": 1e5,
                                "make_peakvels": False},)
