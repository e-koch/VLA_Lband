
'''
For comparison, smooth and regrid to the 14B data.
'''

import os

from cube_analysis.reprojection import reproject_cube

from paths import (data_path, fourteenB_wGBT_HI_file_dict,
                   seventeenB_02kms_wGBT_HI_file_dict,
                   seventeenB_HI_data_02kms_wGBT_path)


out_folder = seventeenB_HI_data_02kms_wGBT_path("14B_match",
                                                no_check=True)

if not os.path.exists(out_folder):
    os.mkdir(out_folder)

out_name = seventeenB_02kms_wGBT_HI_file_dict['Cube'].split("/")[-1].rstrip(".fits") + \
    ".14B_match.fits"

reproject_cube(seventeenB_02kms_wGBT_HI_file_dict['Cube'],
               fourteenB_wGBT_HI_file_dict['Cube'],
               out_name,
               output_folder=out_folder,
               save_spectral=False,
               is_huge=True,
               reproject_type='spatial',
               common_beam=True,
               verbose=True,
               chunk=40)
