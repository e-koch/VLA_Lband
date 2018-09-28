
'''
Plot the HI Zeroth moment with CO contours overlaid.
'''

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from radio_beam import Beam
from astropy.wcs.utils import proj_plane_pixel_area
import seaborn as sb

from paths import (fourteenB_wGBT_HI_file_dict, allfigs_path,
                   iram_co21_14B088_data_path)
from constants import (hi_freq, hi_coldens_Kkms)
from plotting_styles import default_figure, twocolumn_figure
from galaxy_params import gal_feath as gal

default_figure()

cosinc = np.cos(gal.inclination.to(u.rad)).value

moment0 = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]
moment0_wcs = WCS(moment0.header)

radius = gal.radius(header=moment0.header)

beam = Beam.from_fits_header(moment0.header)

# Convert to K km s and correct for disk inclination.
moment0_Kkm_s = beam.jtok(hi_freq).value * (moment0.data / 1000.) * cosinc
moment0_coldens = moment0_Kkm_s * hi_coldens_Kkms.value

pixscale = np.sqrt(proj_plane_pixel_area(moment0_wcs))

# Use the reprojected version
co_moment0 = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits"))[0]

co_noise_map = fits.open(iram_co21_14B088_data_path("m33.rms.14B-088_HI.fits"))[0]

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
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')
# ax.add_patch(beam.ellipse_to_plot(int(0.05 * moment0.shape[0]),
#                                   int(0.05 * moment0.shape[1]), pixscale))

length = (1. * u.kpc / (840 * u.kpc)).to(u.deg, u.dimensionless_angles())
length_pix = length.value / np.abs(moment0.header['CDELT2'])
ax.plot([200, 200], [200 - length_pix / 2., 200 + length_pix / 2.], 'k',
        linewidth=2)
ax.text(80, 200,
        "1 kpc", color='k', va='center')

ax.contour(np.isfinite(co_moment0.data).astype(float), levels=[0.5],
           colors=sb.color_palette('viridis')[:1], alpha=0.45)
ax.contour(co_moment0.data,
           levels=[900, 1400, 1900, 2400],
           cmap='viridis', alpha=0.45)

# CO levels at: 877,  1370,  1862,  2354 K m/s
# <Quantity [  5.8759,  9.179 , 12.4754, 15.7718] solMass / pc2>

cbar = plt.colorbar(im)
cbar.set_label(r"HI Column Density (cm$^{-2}$)")

plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21.pdf"))
plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21.png"))

plt.close()

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
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')
# ax.add_patch(beam.ellipse_to_plot(int(0.05 * moment0.shape[0]),
#                                   int(0.05 * moment0.shape[1]), pixscale))

ax.plot([200, 200], [200 - length_pix / 2., 200 + length_pix / 2.], 'k',
        linewidth=2)
ax.text(80, 200,
        "1 kpc", color='k', va='center')

ax.contour(np.isfinite(co_moment0.data).astype(float), levels=[0.5],
           colors=sb.color_palette('viridis')[:1], alpha=0.45)
ax.contour(co_moment0.data,
           levels=[900, 1400, 1900, 2400],
           cmap='viridis', alpha=0.45)
ax.contour(np.isfinite(co_noise_map.data),
           levels=[0.5], colors=sb.color_palette()[::-1],
           linewidths=[3])

cbar = plt.colorbar(im)
cbar.set_label(r"HI Column Density (cm$^{-2}$)")

plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21_map_edge.pdf"))
plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21_map_edge.png"))

plt.close()

# And with the 7 kpc cut-off region shown

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
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')
# ax.add_patch(beam.ellipse_to_plot(int(0.05 * moment0.shape[0]),
#                                   int(0.05 * moment0.shape[1]), pixscale))

ax.plot([200, 200], [200 - length_pix / 2., 200 + length_pix / 2.], 'k',
        linewidth=2)
ax.text(80, 200,
        "1 kpc", color='k', va='center')

ax.contour(radius < 7 * u.kpc, levels=[0.5], colors=[sb.color_palette()[4]],
           linewidths=[4], alpha=0.5, linestyles='--')

ax.contour(np.isfinite(co_moment0.data).astype(float), levels=[0.5],
           colors=sb.color_palette('viridis')[:1], alpha=0.45)
ax.contour(co_moment0.data,
           levels=[900, 1400, 1900, 2400],
           cmap='viridis', alpha=0.45)
ax.contour(np.isfinite(co_noise_map.data),
           levels=[0.5], colors=sb.color_palette()[::-1],
           linewidths=[3])

cbar = plt.colorbar(im)
cbar.set_label(r"HI Column Density (cm$^{-2}$)")

plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21_map_edge_7kpc.pdf"))
plt.savefig(allfigs_path("HI_maps/coldens_map_14B088_w_CO21_map_edge_7kpc.png"))

plt.close()

default_figure()
