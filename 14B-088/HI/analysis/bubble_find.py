
import os
from spectral_cube import SpectralCube
import astropy.units as u

from basics.bubble_segment3D import BubbleFinder

'''
Create the bubble catalogue of M33 with the 14B-088 map
'''

data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"

cube = \
    SpectralCube.read(os.path.join(data_path,
                                    "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))

# The first few channels have some emission in them. It's partially M33,
# partially galactic HI. So estimate the noise level from the last channel
# which is pretty much noise only. It gives 1.8 mJy/bm, which is right on.
bub_find = BubbleFinder(cube, keep_threshold_mask=True,
                        empty_channel=cube.shape[0] - 1)

bub_find.get_bubbles(verbose=True, overlap_frac=0.5, multiprocess=True,
                     refit=False, nsig=2.0, cut_val=0.5, min_channels=5,
                     distance=0.84 * u.Mpc, nprocesses=None)

bub_find.visualize_bubbles()
