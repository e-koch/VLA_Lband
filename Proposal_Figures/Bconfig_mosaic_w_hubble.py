
'''
Plot the proposed VLA mosaics along with Hubble's coverage.
'''

import sys
import os
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from aplpy import FITSFigure

# Add parent directory to import the directory with the data
sys.path.insert(1, os.path.join(sys.path[0], '..'))

# from plotting_styles import twocolumn_figure

# twocolumn_figure()

hi_freq = 1.42040575177 * u.GHz

# Using the zeroth moment from the 'old' just because it has a lower pb cut
# and looks qualitatively the same
old_data_path = "/mnt/MyRAID/M33/VLA/14B-088/HI/full_imaging"
moment0_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits"
moment0 = fits.open(os.path.join(old_data_path, moment0_name))[0]
moment0_wcs = WCS(moment0.header)

fig = FITSFigure(os.path.join(old_data_path, moment0_name))
fig.show_grayscale()
fig.show_regions("m33_bconfig_7point.reg")
fig.show_regions("m33_bconfig_3point.reg")
fig.show_regions("m33_bricks.reg")

fig.hide_axis_labels()
fig.hide_xtick_labels()
fig.hide_ytick_labels()
plt.tight_layout()

output_path = os.path.expanduser("~/Dropbox/Various Plots/")

fig.save(os.path.join(output_path, "M33_HI_14B088_Bconfig_hubble.pdf"))
fig.save(os.path.join(output_path, "M33_HI_14B088_Bconfig_hubble.png"))

fig.close()
