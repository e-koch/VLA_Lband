
'''
Plot the HI Zeroth moment with CO contours overlaid.
'''

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import astropy.units as u
from radio_beam import Beam
from astropy.wcs.utils import proj_plane_pixel_area

from paths import (fourteenB_HI_data_path, paper1_figures_path,
                   iram_co21_data_path)
from constants import (moment0_name, hi_freq, hi_coldens_Kkms)
from plotting_styles import default_figure, twocolumn_figure
from galaxy_params import gal

default_figure()

moment0 = fits.open(fourteenB_HI_data_path(moment0_name))[0]
moment0_wcs = WCS(moment0.header)

beam = Beam.from_fits_header(moment0.header)

moment0_Kkm_s = beam.jtok(hi_freq).value * moment0.data / 1000.
moment0_coldens = moment0_Kkm_s * hi_coldens_Kkms.value

pixscale = np.sqrt(proj_plane_pixel_area(moment0_wcs))

# Use the reprojected version
co_moment0 = fits.open(iram_co21_data_path("m33.ico.hireprojection.fits"))[0]
# Mask above ~6 kpc
co_moment0 = co_moment0.data
co_moment0[gal.radius(header=moment0.header) > 6 * u.kpc] = np.NaN

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

ax.contour(co_moment0,
           levels=np.linspace(2.0, np.nanpercentile(co_moment0, 99.95), 4),
           cmap='viridis')

cbar = plt.colorbar(im)
cbar.set_label(r"Integrated Intensity (K km s$^{-1}$)")

# raw_input("?")

plt.savefig(paper1_figures_path("zeroth_moment_map_14B088_w_CO21.pdf"))
plt.savefig(paper1_figures_path("zeroth_moment_map_14B088_w_CO21.png"))

plt.close()
