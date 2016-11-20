
import os
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as p
import numpy as np

from analysis.paths import fourteenB_HI_data_path, paper1_figures_path

'''
Plot the residual velocities from disk fit, and from subtracting the smoothed version.
'''

direc = fourteenB_HI_data_path("diskfit_noasymm_nowarp_output/", no_check=True)

res = fits.open(os.path.join(direc, "rad.res.fits"))[0]
res_fit = fits.open(os.path.join(direc, "rad.fitres.fits"))[0]

ax2 = p.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(res.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = p.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

p.draw()

p.savefig(paper1_figures_path("M33_residual_velocity_diskfit.pdf"))
p.savefig(paper1_figures_path("M33_residual_velocity_diskfit.png"))

raw_input("Next plot?")
p.clf()

ax2 = p.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = p.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

p.draw()

p.savefig(paper1_figures_path("M33_residual_velocity_diskfit_fit.pdf"))
p.savefig(paper1_figures_path("M33_residual_velocity_diskfit_fit.png"))
