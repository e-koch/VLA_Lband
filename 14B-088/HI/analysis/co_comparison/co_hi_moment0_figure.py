
'''
Plot the HI Zeroth moment with CO contours overlaid.
'''

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from radio_beam import Beam
from astropy.wcs.utils import proj_plane_pixel_area
import seaborn as sb

from paths import (fourteenB_HI_file_dict, allfigs_path,
                   iram_co21_14B088_data_path)
from constants import (hi_freq, hi_coldens_Kkms)
from plotting_styles import default_figure, twocolumn_figure
# from galaxy_params import gal

default_figure()

moment0 = fits.open(fourteenB_HI_file_dict["Moment0"])[0]
moment0_wcs = WCS(moment0.header)

beam = Beam.from_fits_header(moment0.header)

moment0_Kkm_s = beam.jtok(hi_freq).value * moment0.data / 1000.
moment0_coldens = moment0_Kkm_s * hi_coldens_Kkms.value

pixscale = np.sqrt(proj_plane_pixel_area(moment0_wcs))

# Use the reprojected version
co_moment0 = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits"))[0]

twocolumn_figure(fig_ratio=0.95, font_scale=1.2)

ax = plt.subplot(111, projection=moment0_wcs)
im = ax.imshow(moment0_coldens,
               origin='lower',
               interpolation='nearest',
               alpha=0.85,
               norm=ImageNormalize(vmin=-0.001,
                                   vmax=np.nanmax(moment0_coldens),
                                   stretch=AsinhStretch()))
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
ax.add_patch(beam.ellipse_to_plot(int(0.05 * moment0.shape[0]),
                                  int(0.05 * moment0.shape[1]), pixscale))

ax.contour(np.isfinite(co_moment0.data).astype(float), levels=[0.5],
           colors=sb.color_palette('viridis')[:1])
ax.contour(co_moment0.data,
           levels=np.linspace(np.nanpercentile(co_moment0.data, 70),
                              np.nanpercentile(co_moment0.data, 95), 4),
           cmap='viridis')

# CO levels at: 877,  1370,  1862,  2354 K m/s

cbar = plt.colorbar(im)
cbar.set_label(r"HI Column Density (cm$^{-2}$)")

plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21.pdf"))
plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21.png"))

plt.close()

default_figure()
