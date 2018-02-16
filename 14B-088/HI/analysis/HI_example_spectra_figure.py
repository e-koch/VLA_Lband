
'''
Make a figure showing off some of the different spectra.
Compare archival VLA, new VLA, and VLA+GBT.
'''

from os.path import join as osjoin
from spectral_cube import SpectralCube, Projection
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from radio_beam import Beam
import os
import seaborn as sb
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import numpy as np

from paths import (fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   data_path, allfigs_path)
from plotting_styles import twocolumn_figure, default_figure
from constants import hi_freq, hi_coldens_Kkms
from galaxy_params import gal_feath as gal


figure_folder = allfigs_path("HI_maps")
if not os.path.exists(allfigs_path(figure_folder)):
    os.mkdir(allfigs_path(figure_folder))

arch_cube = SpectralCube.read(osjoin(data_path, "VLA/AT0206/old_imaging/m33_hi_matchto_14B088.fits"))
arch_cube_orig = SpectralCube.read(osjoin(data_path, "VLA/AT0206/old_imaging/m33_hi.fits"))
vla_cube = SpectralCube.read(fourteenB_HI_file_dict['Cube'])
comb_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])

peak_vels = fits.open(fourteenB_wGBT_HI_file_dict['PeakVels'])[0]
peak_temps = fits.open(fourteenB_wGBT_HI_file_dict['PeakTemp'])[0]

# example_posns = [(924, 445), (986, 505), (624, 955), (651, 735),
#                  (659, 745), (802, 640)]
# example_labels = ["a", "b", "c", "d", "e", "f"]

example_posns = [(986, 505), (924, 445),
                 (659, 745), (802, 640)]
example_labels = ["a", "b", "c", "d",]


col_pal = sb.color_palette()

twocolumn_figure()

gs = gridspec.GridSpec(4, 4)

# Plot the column density map with locations of the spectra
moment0_feath = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]
moment0_wcs = WCS(moment0_feath.header)

beam = Beam.from_fits_header(moment0_feath.header)
moment0_Kkm_s_feath = beam.jtok(hi_freq).value * moment0_feath.data / 1000.
moment0_coldens_feath = moment0_Kkm_s_feath * hi_coldens_Kkms.value

ax1 = plt.subplot(gs[:, :2], projection=moment0_wcs)
im = ax1.imshow(moment0_coldens_feath,
                origin='lower',
                interpolation='nearest',
                norm=ImageNormalize(vmin=-0.001 * hi_coldens_Kkms.value,
                                    vmax=np.nanmax(moment0_coldens_feath),
                                    stretch=AsinhStretch()))
ax1.set_ylabel("")
ax1.set_xlabel("")

# Add a box for the location of the N plume figure
ax1.add_patch(patches.Rectangle((500, 1280), 250, 250, fill=False,
                                edgecolor=col_pal[2], linewidth=2,
                                linestyle='dashed'))

# Add a couple radial contours
rad = gal.radius(header=moment0_feath.header).to(u.kpc)
ax1.contour(rad.value, levels=[2, 4], colors=col_pal[:2], linestyle='dotted',
            linewidth=2)

col_pal = sb.color_palette()
lab_posn = [(100, 1100), (100, 800), (1130, 600), (1100, 1000)]
props = dict(boxstyle='round', alpha=1.0, facecolor='w', pad=0.3)
for j, (y, x) in enumerate(example_posns):
    # ax1.scatter(x, y, color=col_pal[j], edgecolor='w',
    #             marker='D', s=70)
    ax1.annotate("({})".format(example_labels[j]), xy=(x, y),
                 xytext=lab_posn[j],
                 fontsize=14, color=col_pal[j],
                 bbox=props,
                 arrowprops=dict(facecolor=col_pal[j], shrink=0.01,
                                 width=4),
                 verticalalignment='center',
                 horizontalalignment='center')

ax2 = plt.subplot(gs[0, 2:])
ax3 = plt.subplot(gs[1, 2:])
ax4 = plt.subplot(gs[2, 2:])
ax5 = plt.subplot(gs[3, 2:])

for i, ax in enumerate([ax2, ax3, ax4, ax5]):

    y, x = example_posns[i]
    print(example_posns[i])

    # All convolved to same beam
    # Though the archival cube's units are still in Jy/bm for the original
    # beam size
    arch_spec = arch_cube[:, y, x].value * arch_cube_orig.beam.jtok(hi_freq).value * u.K
    vla_spec = vla_cube[:, y, x].value * arch_cube.beam.jtok(hi_freq).value * u.K
    comb_spec = comb_cube[:, y, x].value * arch_cube.beam.jtok(hi_freq).value * u.K

    ax.plot((arch_cube.spectral_axis.value - peak_vels.data[y, x]) / 1000.,
            arch_spec.value, color='gray', alpha=0.8, label='Archival',
            drawstyle='steps-mid')
    # ax.plot(vla_cube.spectral_axis.value / 1000., vla_spec, label='VLA',
    #         drawstyle='steps-mid')
    ax.plot((comb_cube.spectral_axis.value - peak_vels.data[y, x]) / 1000.,
            comb_spec, label='VLA+GBT', color='k',
            drawstyle='steps-mid')
    ax.grid()
    ax.axhline(0, color='k', alpha=0.5, linestyle='--')

    low_lim = max(vla_cube.spectral_axis.min().value,
                  peak_vels.data[y, x] - 50000)
    high_lim = min(vla_cube.spectral_axis.max().value,
                   peak_vels.data[y, x] + 50000)
    ax.set_xlim([-50000 / 1000., 50000 / 1000.])

    # if i == 0:
    #     ax.legend(loc='upper left', frameon=True)

    # if i == 1 or i == 3 or i == 5:
    # if i == 2 or i == 3:
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')

    # if i == 2:
    # if i == 2 or i == 3:
    ax.set_ylabel(r"$T_{\rm HI}$ (K)")

    # if i == 4 or i == 5:
    if i == 3:
        ax.set_xlabel(r"Velocity (km/s)")
    else:
        ax.axes.xaxis.set_ticklabels([])

    ax.text(0.9, 0.8, "({})".format(example_labels[i]),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            bbox={"boxstyle": "square", "facecolor": "w"},
            fontsize=14, weight='bold',
            # Match the color with the symbols in the column density figure
            color=col_pal[i])

plt.tight_layout()

plt.savefig(osjoin(figure_folder, "example_HI_spectra.png"))
plt.savefig(osjoin(figure_folder, "example_HI_spectra.pdf"))

plt.close()

# peak_temps_proj = Projection.from_hdu(peak_temps)
# peak_temps_proj.quicklook()
# dec, ra = comb_cube.spatial_coordinate_map
# for j, (y, x) in enumerate(example_posns):
#     peak_temps_proj.FITSFigure._ax1.scatter(x, y, color=col_pal[j], marker='x')

default_figure()
