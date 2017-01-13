
import os
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as p
import numpy as np

from analysis.paths import fourteenB_HI_data_path, paper1_figures_path
from analysis.plotting_styles import default_figure, twocolumn_twopanel_figure
from analysis.constants import moment1_name

'''
Plot the residual velocities from disk fit, and from subtracting the
smoothed version.
'''

direc = fourteenB_HI_data_path("diskfit_noasymm_noradial_nowarp_output/",
                               no_check=True)

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

# raw_input("Next plot?")
p.clf()

ax2 = p.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = p.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

p.draw()

p.savefig(paper1_figures_path("M33_residual_velocity_diskfit_fit.pdf"))
p.savefig(paper1_figures_path("M33_residual_velocity_diskfit_fit.png"))

p.close()

# Make a side-by-side centroid and residual plot.

twocolumn_twopanel_figure(fig_ratio=0.5)

ax2 = p.subplot(122, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)
ax2.set_xlabel("RA (J2000)")
ax2.coords[1].set_ticklabel_visible(False)

cb2 = p.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

moment1 = fits.open(fourteenB_HI_data_path(moment1_name))[0]

ax = p.subplot(121, projection=WCS(moment1.header))
im = ax.imshow(moment1.data / 1000.,
               origin='lower', vmin=-300, vmax=-70,
               interpolation='nearest', cmap='viridis')

ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")

cbar = p.colorbar(im)
cbar.set_label(r"Centroid Velocity (km s$^{-1}$)")

p.savefig(paper1_figures_path("M33_centroid_residualfit.pdf"))
p.savefig(paper1_figures_path("M33_centroid_residualfit.png"))

p.close()

default_figure()
