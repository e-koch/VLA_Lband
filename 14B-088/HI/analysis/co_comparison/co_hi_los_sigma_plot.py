
'''
Highlight a couple of different regions to show that the local variation
amongst the aggregate distribution.
'''

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import Projection
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from radio_beam import Beam
from astropy.wcs.utils import proj_plane_pixel_area
import seaborn as sb
from astropy.table import Table
from corner import hist2d
import matplotlib.patches as patches
import os

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path,
                   iram_co21_14B088_data_path)
from constants import (hi_freq, hi_coldens_Kkms)
from plotting_styles import default_figure, onecolumn_Npanel_figure
from galaxy_params import gal_feath as gal

default_figure()

fig_path = allfigs_path("co_vs_hi")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

col_pal = sb.color_palette()

cosinc = np.cos(gal.inclination.to(u.rad)).value

moment0 = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]
moment0_wcs = WCS(moment0.header)

mom0_proj = Projection.from_hdu(moment0)

beam = Beam.from_fits_header(moment0.header)

# Convert to K km s and correct for disk inclination.
moment0_Kkm_s = beam.jtok(hi_freq).value * (moment0.data / 1000.) * cosinc
moment0_coldens = moment0_Kkm_s * hi_coldens_Kkms.value

pixscale = np.sqrt(proj_plane_pixel_area(moment0_wcs))

# Use the reprojected version
co_moment0 = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits"))[0]

co_noise_map = fits.open(iram_co21_14B088_data_path("m33.rms.14B-088_HI.fits"))[0]

onecolumn_Npanel_figure(N=1.5)

img_slice = (slice(902, 1084), slice(400, 650))
offset = [902, 400]

fig = plt.figure(figsize=(4.4, 5.84))

ax = fig.add_axes((0.15, 0.5, 0.8, 0.45), projection=mom0_proj[img_slice].wcs)
im = ax.imshow(moment0_coldens[img_slice],
               origin='lower',
               interpolation='nearest',
               alpha=0.85,
               norm=ImageNormalize(vmin=-0.001,
                                   vmax=np.nanmax(moment0_coldens),
                                   stretch=AsinhStretch()))

lon = ax.coords[0]
lon.set_major_formatter('hh:mm:ss')


ax.contour(np.isfinite(co_moment0.data[img_slice]).astype(float), levels=[0.5],
           colors=sb.color_palette('viridis')[:1],
           alpha=0.5)
ax.contour(co_moment0.data[img_slice],
           levels=np.linspace(np.nanpercentile(co_moment0.data, 70),
                              np.nanpercentile(co_moment0.data, 95), 4),
           cmap='viridis', alpha=0.5)

reg1 = ((111, 31), (143, 62))
reg2 = ((142, 82), (157, 99))
reg3 = ((90, 132), (138, 170))

ax.add_patch(patches.Rectangle(reg1[0], reg1[1][0] - reg1[0][0],
                               reg1[1][1] - reg1[0][1], fill=True,
                               facecolor=col_pal[3], linewidth=1,
                               edgecolor='k'))

ax.add_patch(patches.Rectangle(reg2[0], reg2[1][0] - reg2[0][0],
                               reg2[1][1] - reg2[0][1], fill=True,
                               facecolor=col_pal[2], linewidth=1,
                               edgecolor='k'))

ax.add_patch(patches.Rectangle(reg3[0], reg3[1][0] - reg3[0][0],
                               reg3[1][1] - reg3[0][1], fill=True,
                               facecolor=col_pal[5], linewidth=1,
                               edgecolor='k'))


tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_co_gaussfit_column_densities_perpix.fits"))

# Don't consider the "bad fits" that are probably due to multiple components
good_pts = np.logical_and(~tab['multicomp_flag_HI'],
                          ~tab['multicomp_flag_CO'])
good_pts = np.logical_and(good_pts,
                          tab["sigma_HI"] > 3800)
# Minimum CO line width of one channel.
good_pts = np.logical_and(good_pts,
                          tab["sigma_CO"] >= 2600)

ax2 = fig.add_axes((0.15, 0.1, 0.8, 0.35))

hist2d(tab['sigma_HI'][good_pts] / 1000.,
       np.array(tab['sigma_CO'])[good_pts] / 1000., bins=13,
       data_kwargs={"alpha": 0.5}, ax=ax2)
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")

# From line_ratios_analysis.py
slope_ratio = 0.56318
plt.plot([4, 12], [4. * slope_ratio, 12. * slope_ratio],
         '--', color=sb.color_palette()[1], linewidth=3,
         label='Ratio Fit', alpha=0.8)

# plt.plot([4, 12], [4, 12], '-.', linewidth=3, alpha=0.6,
#          color=sb.color_palette()[2],
#          label=r'$\sigma_{\rm CO} = \sigma_{\rm HI}$')

ax2.set_ylim([0.5, 8])
ax2.legend(frameon=True, loc='lower right')
ax2.grid()

ax2.set_xlabel(r"$\sigma_{\rm HI}$ (km/s)")
ax2.set_ylabel(r"$\sigma_{\rm CO}$ (km/s)")

# Overplot the points in the three regions


def return_pts_mask(reg, tab, offset):

    x_sel = np.logical_and(tab['xpts'] >= reg[0][0] + offset[1],
                           tab['xpts'] <= reg[1][0] + offset[1])
    y_sel = np.logical_and(tab['ypts'] >= reg[0][1] + offset[0],
                           tab['ypts'] <= reg[1][1] + offset[0])

    return np.logical_and(x_sel, y_sel)


reg1_mask = np.logical_and(good_pts, return_pts_mask(reg1, tab, offset))
reg2_mask = np.logical_and(good_pts, return_pts_mask(reg2, tab, offset))
reg3_mask = np.logical_and(good_pts, return_pts_mask(reg3, tab, offset))

ax2.scatter(tab['sigma_HI'][reg1_mask] / 1000.,
            tab['sigma_CO'][reg1_mask] / 1000., c=col_pal[3],
            marker='D')
ax2.scatter(tab['sigma_HI'][reg2_mask] / 1000.,
            tab['sigma_CO'][reg2_mask] / 1000., c=col_pal[2],
            marker='o')
ax2.scatter(tab['sigma_HI'][reg3_mask] / 1000.,
            tab['sigma_CO'][reg3_mask] / 1000., c=col_pal[5],
            marker='s')

plt.savefig(os.path.join(fig_path, "sigma_HI_vs_H2_w_example_regions.png"))
plt.savefig(os.path.join(fig_path, "sigma_HI_vs_H2_w_example_regions.pdf"))
plt.close()
