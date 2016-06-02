
from aplpy import FITSFigure
import os

mom0_file = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits"

fig = FITSFigure(mom0_file)
fig.show_grayscale()
fig.show_regions(os.path.expanduser("~/Dropbox/code_development/VLA_Lband/14B-088/pointings.reg"))
