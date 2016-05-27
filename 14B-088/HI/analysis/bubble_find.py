
import os
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np

from basics.bubble_segment3D import BubbleFinder
from basics.utils import sig_clip

'''
Create the bubble catalogue of M33 with the 14B-088 map
'''

data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"

cube = \
    SpectralCube.read(os.path.join(data_path,
                                    "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))

sigma = sig_clip(cube[-1].value, nsig=10)

# Create a cruddy linewidth map. This will need to be improved, but that can
# be done with just the saved bubble objects
lwidth = cube.with_mask(cube > 2 * sigma * u.Jy).linewidth_sigma()

galaxy_props = {"center_coord":
                SkyCoord(23.461667, 30.660194,
                         unit=(u.deg, u.deg), frame='fk5'),
                "inclination": 56. * u.deg,
                "position_angle": 201. * u.deg,
                "scale_height": 100. * u.pc}


scales = 3. * np.arange(1, 15, np.sqrt(2))

# The first few channels have some emission in them. It's partially M33,
# partially galactic HI. So estimate the noise level from the last channel
# which is pretty much noise only. It gives 1.8 mJy/bm, which is right on.
bub_find = BubbleFinder(cube, keep_threshold_mask=True,
                        # empty_channel=cube.shape[0] - 1,  # Given above
                        sigma=sigma,
                        galaxy_props=galaxy_props,
                        distance=0.84 * u.Mpc)

# Requiring minimum of 15 channels, similar to the 3 channel requirement with
# THINGS (15 * 0.2 km/s = 3 km/s; 3 * 1.2 km/s = 3.6 km/s)
bub_find.get_bubbles(verbose=True, overlap_frac=0.5, multiprocess=True,
                     refit=False, nsig=1.5, min_corr=0.7, min_overlap=0.8,
                     min_channels=15, nprocesses=None, scales=scales,
                     cube_linewidth=lwidth)

# Save all of the bubbles
folder = os.path.expanduser("~/MyRAID/M33/14B-088/HI/full_imaging/bubbles/")
try:
    os.mkdir(folder)
except OSError:
    pass
bub_find.save_bubbles(folder=folder, name="M33_14B-088")

catalog = bub_find.to_catalog()

catalog.write_table(os.path.join(folder, "bubble_catalog.ecsv"))

# Show the outlines of the bubbles
bub_find.visualize_bubbles()
