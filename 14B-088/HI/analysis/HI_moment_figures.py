
'''
Make figures of the zeroth and first moments.
'''

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from radio_beam import Beam
from astropy.wcs.utils import proj_plane_pixel_area

from paths import fourteenB_HI_data_path, paper1_figures_path
from constants import moment0_name, moment1_name, hi_freq
from plotting_styles import default_figure, twocolumn_figure

default_figure()

moment0 = fits.open(fourteenB_HI_data_path(moment0_name))[0]
moment0_wcs = WCS(moment0.header)

beam = Beam.from_fits_header(moment0.header)

moment0_Kkm_s = beam.jtok(hi_freq).value * moment0.data / 1000.

pixscale = np.sqrt(proj_plane_pixel_area(moment0_wcs))

twocolumn_figure(fig_ratio=0.95, font_scale=1.2)

ax = plt.subplot(111, projection=moment0_wcs)
im = ax.imshow(moment0_Kkm_s,
               origin='lower',
               interpolation='nearest',
               norm=ImageNormalize(vmin=-0.001,
                                   vmax=np.nanmax(moment0_Kkm_s),
                                   stretch=AsinhStretch()))
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
ax.add_patch(beam.ellipse_to_plot(int(0.05 * moment0.shape[0]),
                                  int(0.05 * moment0.shape[1]), pixscale))

cbar = plt.colorbar(im)
cbar.set_label(r"Intensity (K km s$^{-1}$)")

# raw_input("?")

plt.savefig(paper1_figures_path("zeroth_moment_map_14B088.pdf"))
plt.savefig(paper1_figures_path("zeroth_moment_map_14B088.png"))

plt.close()

default_figure()

moment1 = fits.open(fourteenB_HI_data_path(moment1_name))[0]

ax = plt.subplot(111, projection=WCS(moment1.header))
im = ax.imshow(moment1.data / 1000.,
               origin='lower', vmin=-300, vmax=-70,
               interpolation='nearest', cmap='viridis')

ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")

cbar = plt.colorbar(im)
cbar.set_label(r"Centroid Velocity (km s$^{-1}$)")

plt.savefig(paper1_figures_path("centroid_map_14B088.pdf"))
plt.savefig(paper1_figures_path("centroid_map_14B088.png"))

plt.close()
