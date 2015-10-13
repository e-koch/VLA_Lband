
# Remove continuum outliers around M33
# Coordinates can be found in text file

from casa_tools import subtract_outliers

import glob


# Grab the list of mask images
masks = glob.glob("masks/*.mask")
masks.sort()
assert len(masks) == 9

# Arranged roughly in decreasing brightness
outlier_coords = ['J2000 1h32m22.618 +30d44m08.711',
                  'J2000  1h34m39.033 +30d03m50.412',
                  "J2000 1h36m54.049 +30d14m35.373",
                  "J2000 1h34m29.570 +31d03m12.081",
                  "J2000 1h35m13.195 31d10m33.247",
                  "J2000 1h35m16.252 30d54m51.118",
                  "J2000 1h35m08.286 31d05m35.170",
                  "J2000 1h36m17.153 31d24m34.212",
                  "J2000 1h37m08.449 31d22m46.636"]

# Now run the subtraction function
subtract_outliers('14B-088_continuum_I.ms', outlier_coords,
                  threshold='0.1mJy/beam', datacolumn='DATA',
                  field='M33_7_center', interactive=False, cleanup=True,
                  save_space=True, masks=masks)
