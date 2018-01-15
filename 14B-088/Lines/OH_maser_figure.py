
'''
Research note figure on the OH(1665) detection
'''

from spectral_cube import SpectralCube, Projection
from aplpy import FITSFigure
from astropy.io import fits
import pyregion
from os.path import join as osjoin
from astropy.wcs.utils import proj_plane_pixel_area
import astropy.units as u
import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt
import os
import seaborn as sb

from constants import hi_freq
from plotting_styles import onecolumn_figure
from paths import allfigs_path, fourteenB_wGBT_HI_file_dict, data_path, c_path

onecolumn_figure()

cpal = sb.color_palette()

fig_folder = os.path.join(allfigs_path("OH_maser"))
if not os.path.exists(fig_folder):
    os.mkdir(fig_folder)

path = "/mnt/bigdata/ekoch/VLA/14B-088/Lines/OH/OH1665/detection_imaging_1point5km_s"

# Open the OH 1665 narrow line width cube

cube = SpectralCube.read(osjoin(path, "OH1665_14B-088_uniform.image.pbcor.fits"))

# Convolve to a common beam
com_beam = cube.beams.common_beam(epsilon=7e-4)
cube = cube.convolve_to(com_beam)

reg = pyregion.open("/mnt/bigdata/ekoch/VLA/14B-088/Lines/OH/oh1665_maser_uniform.reg")

reg_cube = cube.subcube_from_ds9region(reg)
mom0 = reg_cube.spectral_slab(-220 * u.km / u.s, -200 * u.km / u.s).moment0()

# Integrated intensity over region

num_pix = np.sum(np.isfinite(mom0)) * u.pix
pix_size = proj_plane_pixel_area(mom0.wcs) * u.Unit(mom0.header['CUNIT2'])**2

pix_per_beam = (com_beam.sr.to(u.deg**2) / pix_size) * (u.pix / u.beam)

# Make a total spectrum based on the peak location
sum_spec = reg_cube[:, 3, 4] * u.beam

mad_std = 0.0018596 * u.Jy

onecolumn_figure()

plt.plot(sum_spec.spectral_axis.to(u.km / u.s), sum_spec.value,
         drawstyle='steps-mid')
plt.ylabel("Flux (Jy)", fontsize=13)
plt.xlabel("Velocity (km/s)", fontsize=13)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.xlim(-212 - 12, -212 + 12)
plt.axhline(mad_std.value, color=sb.color_palette()[1], linestyle='--')
plt.tight_layout()
plt.savefig(os.path.join(fig_folder, "OH1665_spectrum.png"))
plt.savefig(os.path.join(fig_folder, "OH1665_spectrum.pdf"))
plt.close()

# Next load up the HI mom0 image
hi_hdu = fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0]

# BUNIT says Jy m/s, but it's really Jy/bm m/s
hi_proj = Projection.from_hdu(hi_hdu)

# Convert to K km/s
hi_Kkms = (hi_proj * (hi_proj.beam.jtok(hi_freq) / u.Jy)).to(u.K * u.km / u.s)

fig = plt.figure(figsize=(4, 4))

hi_fig = FITSFigure(hi_Kkms.hdu, figure=fig)
hi_fig.show_grayscale(invert=True)
# hi_fig.show_colorbar()
# hi_fig.colorbar.set_axis_label_text('Integrated Intensity (K km/s)')
# hi_fig.tick_labels.set_xformat('h:mm:ss')
# hi_fig.tick_labels.set_yformat('dd:mm:ss')
hi_fig.hide_axis_labels()
hi_fig.hide_tick_labels()

x_world = 23.5009
y_world = 30.680
hi_fig.show_markers(x_world, y_world, marker="D",
                    facecolor=cpal[2],
                    s=30, edgecolor=cpal[2],
                    linewidth=3)

plt.tight_layout()

hi_fig.savefig(os.path.join(fig_folder, 'OHmaser_HImap.png'))
hi_fig.savefig(os.path.join(fig_folder, 'OHmaser_HImap.pdf'))
hi_fig.close()

# Now a close-up in the Halpha with labeled regions.

halpha_hdu = fits.open(os.path.join(data_path, "Halpha/ha.fits"))

halpha_proj = Projection.from_hdu(halpha_hdu)

regions_slice = (slice(8120, 8280), slice(4010, 4170))

fig = plt.figure(figsize=(5, 4))

halp_fig = FITSFigure(halpha_proj[regions_slice].hdu, figure=fig)
halp_fig.show_grayscale(invert=True, stretch='arcsinh')
# halp_fig.show_colorbar()
# halp_fig.colorbar.set_axis_label_text('UNIT (K km/s)')
halp_fig.tick_labels.set_xformat('h:mm:ss')
halp_fig.tick_labels.set_yformat('dd:mm:ss')
halp_fig.set_tick_labels_size(12)
halp_fig.hide_axis_labels()
halp_fig.show_regions(os.path.join(c_path, 'Lines/ohbeam_with_hiiregions.reg'))

# I want some nicer colors. Hack into the regions set and change to seaborn
# colors
regions = halp_fig.get_layer('region_set_1')
regions.artistlist[0].set_edgecolor(cpal[0])
regions.artistlist[1].set_edgecolor(cpal[0])
regions.artistlist[2].set_edgecolor(cpal[1])
regions.artistlist[3].set_edgecolor(cpal[1])
regions.artistlist[4].set_edgecolor(cpal[1])
regions.artistlist[5].set_edgecolor(cpal[2])

halp_fig.show_arrows(23.505, 30.67995, -0.004, 0, color=cpal[0])
halp_fig.show_arrows(23.5009, 30.6769, 0, 0.002, color=cpal[0])

halp_fig.add_label(23.505, 30.67995, "C1-1", color=cpal[0], fontsize=13,
                   weight='bold',
                   bbox={"boxstyle": "round", "facecolor": "w"})
halp_fig.add_label(23.5009, 30.6769, "C1-2", color=cpal[0], fontsize=13,
                   weight='bold',
                   bbox={"boxstyle": "round", "facecolor": "w"})

halp_fig.savefig(os.path.join(fig_folder, 'OHmaser_Halpmap_wlabels.png'))
halp_fig.savefig(os.path.join(fig_folder, 'OHmaser_Halpmap_wlabels.pdf'))
halp_fig.close()
