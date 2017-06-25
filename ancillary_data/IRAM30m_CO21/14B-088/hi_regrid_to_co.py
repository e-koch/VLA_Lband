
from cube_analysis.reprojection import reproject_cube

from paths import (iram_co21_data_path, fourteenB_HI_data_path,
                   fourteenB_HI_data_wGBT_path, fourteenB_HI_file_dict,
                   fourteenB_wGBT_HI_file_dict)

'''
Regrid the 14B-088 HI data to match the CO in the celestial dimensions.
'''

reproject_cube(fourteenB_HI_file_dict["Cube"],
               iram_co21_data_path("m33.co21_iram.fits"),
               "M33_14B-088_HI.clean.image.co21_iram_regrid.fits",
               output_folder=fourteenB_HI_data_path("", no_check=True),
               reproject_type='spatial',
               save_spectral=False,
               is_huge=True,
               common_beam=False,
               verbose=True,
               chunk=80)


reproject_cube(fourteenB_wGBT_HI_file_dict["Cube"],
               iram_co21_data_path("m33.co21_iram.fits"),
               "M33_14B-088_HI.clean.image.GBT_feathered.co21_iram_regrid.fits",
               output_folder=fourteenB_HI_data_wGBT_path("", no_check=True),
               reproject_type='spatial',
               save_spectral=False,
               is_huge=True,
               common_beam=False,
               verbose=True,
               chunk=80)

# Make version reprojected in spectral dim as well.

reproject_cube(fourteenB_HI_file_dict["Cube"],
               iram_co21_data_path("m33.co21_iram.fits"),
               "M33_14B-088_HI.clean.image.co21_iram_regrid_full.fits",
               output_folder=fourteenB_HI_data_path("", no_check=True),
               reproject_type='all',
               save_spectral=False,
               is_huge=True,
               common_beam=False,
               verbose=True,
               chunk=80)


reproject_cube(fourteenB_wGBT_HI_file_dict["Cube"],
               iram_co21_data_path("m33.co21_iram.fits"),
               "M33_14B-088_HI.clean.image.GBT_feathered.co21_iram_regrid_full.fits",
               output_folder=fourteenB_HI_data_wGBT_path("", no_check=True),
               reproject_type='all',
               save_spectral=False,
               is_huge=True,
               common_beam=False,
               verbose=True,
               chunk=80)
