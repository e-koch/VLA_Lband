
'''
Make a figure showing off some of the different spectra.
Compare archival VLA, new VLA, and VLA+GBT.
'''

from os.path import join as osjoin
from spectral_cube import SpectralCube, Projection
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
import os
from os.path import join as osjoin
import seaborn as sb

from paths import (fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   data_path, allfigs_path)
from plotting_styles import twocolumn_figure, default_figure
from constants import hi_freq


figure_folder = allfigs_path("HI_maps")
if not os.path.exists(allfigs_path(figure_folder)):
    os.mkdir(allfigs_path(figure_folder))

arch_cube = SpectralCube.read(osjoin(data_path, "VLA/AT0206/old_imaging/m33_hi_matchto_14B088.fits"))
arch_cube_orig = SpectralCube.read(osjoin(data_path, "VLA/AT0206/old_imaging/m33_hi.fits"))
vla_cube = SpectralCube.read(fourteenB_HI_file_dict['Cube'])
comb_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])

peak_vels = fits.open(fourteenB_wGBT_HI_file_dict['PeakVels'])[0]
peak_temps = fits.open(fourteenB_wGBT_HI_file_dict['PeakTemp'])[0]

example_posns = [(924, 445), (986, 505), (624, 955), (651, 735),
                 (659, 745), (802, 640)]
example_labels = ["a", "b", "c", "d", "e", "f"]

col_pal = sb.color_palette()

twocolumn_figure()

fig, axs = plt.subplots(3, 2)

for i, ax in enumerate(axs.ravel()):

    y, x = example_posns[i]
    print(example_posns[i])

    # All convolved to same beam
    # Though the archival cube's units are still in Jy/bm for the original
    # beam size
    arch_spec = arch_cube[:, y, x].value * arch_cube_orig.beam.jtok(hi_freq).value * u.K
    vla_spec = vla_cube[:, y, x].value * arch_cube.beam.jtok(hi_freq).value * u.K
    comb_spec = comb_cube[:, y, x].value * arch_cube.beam.jtok(hi_freq).value * u.K

    ax.plot(arch_cube.spectral_axis.value / 1000., arch_spec.value, label='Archival',
            drawstyle='steps-mid')
    ax.plot(vla_cube.spectral_axis.value / 1000., vla_spec, label='VLA',
            drawstyle='steps-mid')
    ax.plot(comb_cube.spectral_axis.value / 1000., comb_spec, label='VLA+GBT',
            drawstyle='steps-mid')
    ax.grid()
    ax.axhline(0, color='k', alpha=0.5, linestyle='--')
    if i == 0:
        ax.legend(loc='upper left', frameon=True)

    if i == 1 or i == 3 or i == 5:
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()
        ax.yaxis.set_ticks_position('both')

    if i == 2:
        ax.set_ylabel(r"$T_{\rm HI}$ (K)")

    if i == 4 or i == 5:
        ax.set_xlabel(r"Velocity (km/s)")

    low_lim = max(vla_cube.spectral_axis.min().value,
                  peak_vels.data[y, x] - 50000)
    high_lim = min(vla_cube.spectral_axis.max().value,
                   peak_vels.data[y, x] + 50000)
    ax.set_xlim([low_lim / 1000., high_lim / 1000.])

    ax.text(0.9, 0.8, "({})".format(example_labels[i]),
            horizontalalignment='center',
            verticalalignment='center',
            transform=ax.transAxes,
            bbox={"boxstyle": "square", "facecolor": "w"},
            fontsize=14, weight='bold',
            # Match the color with the symbols in the column density figure
            color=col_pal[i])

plt.tight_layout()

fig.savefig(osjoin(figure_folder, "example_HI_spectra.png"))
fig.savefig(osjoin(figure_folder, "example_HI_spectra.pdf"))

plt.close()

# peak_temps_proj = Projection.from_hdu(peak_temps)
# peak_temps_proj.quicklook()
# dec, ra = comb_cube.spatial_coordinate_map
# for j, (y, x) in enumerate(example_posns):
#     peak_temps_proj.FITSFigure._ax1.scatter(x, y, color=col_pal[j], marker='x')

default_figure()