
'''
Create a version of Andreas' cube to match the spatial resolution of the
14B-088 data.
'''


from os.path import join as osjoin

from paths import (fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   data_path)

from cube_analysis.reprojection import reproject_cube

reproject_cube(osjoin(data_path, "VLA/AT0206/old_imaging/m33_hi.fits"),
               fourteenB_HI_file_dict['Cube'],
               osjoin(data_path, "VLA/AT0206/old_imaging/m33_hi_matchto_14B088.fits"),
               reproject_type='spatial',
               common_beam=True)
