
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import numpy as np
import matplotlib.pyplot as p
import os
from astropy.io import fits

'''
Investigating skewness and kurtosis in the 14B-088 cube.
'''

data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"

cube = \
    SpectralCube.read(os.path.join(data_path,
                                    "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))

mask = fits.getdata(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked_source_mask.fits"))
mask = mask.astype(bool)

cube = cube.with_mask(mask)

# Skewness: Converting from m^3/s^3 to km/s
mom3 = cube.moment(order=3) ** (1/3.) / 1000.

# mom3_flattened = \
#     Projection(np.arctan(mom3.value / np.nanpercentile(mom3.value, 85)),
#                wcs=mom3.wcs)

# Ranges from ~-5e13 to 1e13. Visualize by applying an arctan transform
# p.imshow(mom3_flattened.value,
#          origin='lower', cmap='RdGy')

# Save the whole thing, a close-up of 604 and the NArm, the Northern "plume",
# and the Southern Arm
# mom3_flattened.quicklook()
mom3.quicklook()

# Kurtosis: Subtract 3 to center on 0 (assuming a Gaussian)
# Then take 4th root to put in units of km/s
mom4 = (cube.moment(order=4) - 3) ** 0.25 / 1000.

# Interesting regions:
# Northern HI infall: [1373, 630], [1373, 640], [1373, 650]
# progression across the flip in the peaks

# SW of galactic centre: [712 580]

# Southern Arm

# NGC 604 [980, 455], [970, 455], [960, 455], [950, 455], [940, 455],
# [930, 455], [920, 455]
