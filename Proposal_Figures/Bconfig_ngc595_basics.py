
'''
Highlight the HI bubble, corresponding to the Giant HII region NGC 595.
Also show the Kitt Peak Halpha map.

Based on the combined 14B-088 & AT0206 cube.
'''

import os
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.io import fits
import astropy.wcs as wcs
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from reproject import reproject_interp

from basics import BubbleFinder
from basics.utils import sig_clip

data_path = "/mnt/MyRAID/M33/VLA/14B-088/combined_HI/"

cube = \
    SpectralCube.read(os.path.join(data_path,
                                   "M33_14B-088_AT0206_HI.clean.image.fits"))

# Lots of crap on the edge since the B-config fields don't cover all of the C
pbcov = fits.open(os.path.join(data_path,
                               "M33_HI_14B-088_AT0206_pbcov.fits"))[0].data[0]
cube = cube.with_mask(pbcov > 0.6).minimal_subcube()

# This channel is nearly empty of signal.
sigma = sig_clip(cube[-1].value, nsig=10)

# Create a cruddy linewidth map. This will need to be improved, but that can
# be done with just the saved bubble objects
lwidth = cube.with_mask(cube > 2 * sigma * u.Jy).linewidth_sigma()

# Lazy zeroth moment
mom0 = cube.with_mask(cube > 2 * sigma * u.Jy).moment0()

galaxy_props = {"center_coord":
                SkyCoord(23.461667, 30.660194,
                         unit=(u.deg, u.deg), frame='fk5'),
                "inclination": 56. * u.deg,
                "position_angle": 201. * u.deg,
                "scale_height": 100. * u.pc}


# Let's just limit to small scales for sake of the visualization.
scales = 4. * np.arange(1, 15, np.sqrt(2))[:5]

# The first few channels have some emission in them. It's partially M33,
# partially galactic HI. So estimate the noise level from the last channel
# which is pretty much noise only. It gives 1.8 mJy/bm, which is right on.
bub_find = BubbleFinder(cube, keep_threshold_mask=True,
                        sigma=sigma,
                        galaxy_props=galaxy_props,
                        distance=0.84 * u.Mpc)

# Requiring minimum of 15 channels, similar to the 3 channel requirement with
# LITTLE THINGS (15 * 0.2 km/s = 3 km/s; 3 * 1.2 km/s = 3.6 km/s)
bub_find.get_bubbles(verbose=True, overlap_frac=0.5, multiprocess=False,
                     refit=False, nsig=1.5, min_corr=0.6, min_overlap=0.6,
                     global_corr=0.3, min_channels=3, nprocesses=None,
                     scales=scales, cube_linewidth=lwidth,
                     save_regions=False)



# Make a pretty plot. 2 subplots of HI and Halpha

cmap = cm.binary
cmap.set_bad('w', 1.)

halpha = fits.open("/home/eric/MyRAID/M33/ha/ha.fits")[0]
co_int = fits.open("/home/eric/MyRAID/M33/co21/m33.ico.fits")[0]

# Regrid to the HI.
co_wcs = wcs.WCS(co_int.header).dropaxis(-1).dropaxis(-1)
co_regrid = reproject_interp((co_int.data.squeeze(), co_wcs), mom0.header)[0]

halpha_regrid = reproject_interp(halpha, mom0.header)[0]

ngc595_slice = [[1530, 1620], [1200, 1290]]

ax1 = plt.subplot(121)
scaled_hi = np.arctan(mom0.value / np.nanpercentile(mom0.value, 90))
bub_find.visualize_bubbles(moment0=scaled_hi, ax=ax1)
ax1.imshow(scaled_hi, cmap=cmap)
ax1.set_xlim(ngc595_slice[1])
ax1.set_ylim(ngc595_slice[0])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.text(1210, 1610, "HI", fontsize=30,
         bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

ax2 = plt.subplot(122)
scaled_halpha = \
    np.arctan(halpha_regrid / np.nanpercentile(halpha_regrid, 99))
bub_find.visualize_bubbles(moment0=scaled_halpha, ax=ax2)
ax2.imshow(scaled_halpha, cmap=cmap, origin='lower')

ax2.contour(co_regrid, colors='g', levels=[2.0, 4.0, 5.0])
ax2.set_xlim(ngc595_slice[1])
ax2.set_ylim(ngc595_slice[0])
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.text(1210, 1610, r"H$\alpha$ with CO(2-1) contours", fontsize=30,
         bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

plt.tight_layout()

plt.savefig(os.path.expanduser("~/Dropbox/Various Plots/Proposals/BaSiCs_ngc595_at0206_14b088_combine.pdf"))
plt.savefig(os.path.expanduser("~/Dropbox/Various Plots/Proposals/BaSiCs_ngc595_at0206_14b088_combine.png"))
