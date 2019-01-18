
'''
Create a cube that is spatially-matched and convolved to the IRAM 30-m
CO(2-1) cube.
'''

import os
from os.path import join as osjoin

from cube_analysis.reprojection import reproject_cube
from cube_analysis.run_pipe import run_pipeline

from paths import (seventeenB_HI_data_1kms_wGBT_path,
                   seventeenB_1kms_wGBT_HI_file_dict, iram_co21_data_path)


out_folder = seventeenB_HI_data_1kms_wGBT_path("iram_co21_match",
                                               no_check=True)

if not os.path.exists(out_folder):
    os.mkdir(out_folder)

out_name = seventeenB_1kms_wGBT_HI_file_dict['Cube'].split("/")[-1].rstrip(".fits") + \
    ".iram30m_spatialmatch.fits"

reproject_cube(seventeenB_1kms_wGBT_HI_file_dict['Cube'],
               iram_co21_data_path("m33.co21_iram.fits"),
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
