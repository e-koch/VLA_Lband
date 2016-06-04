
import os
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as p
import numpy as np

'''
Plot the residual velocities from disk fit, and from subtracting the smoothed version.
'''

direc = "/home/eric/MyRAID/M33/14B-088/HI/full_imaging/diskfit_noasymm_nowarp_output/"

# res = fits.open(os.path.join(direc, "rad.res.fits"))[0]
res_fit = fits.open(os.path.join(direc, "rad.fitres.fits"))[0]

data = res_fit.data

# ax1 = p.subplot(121, projection=WCS(res.header))
# im1 = ax1.imshow(res.data / 1000., origin='lower', cmap=p.cm.gray_r)

# cb = p.colorbar(im1)
# cb.set_label("km s$^{-1}$")

ax2 = p.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = p.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

p.draw()
