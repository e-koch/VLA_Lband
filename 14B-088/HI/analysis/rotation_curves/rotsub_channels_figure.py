import numpy as np
import os
import matplotlib.pyplot as p
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import average_beams
from astropy import coordinates
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy import log
import aplpy
from astropy.visualization import SqrtStretch, AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import warnings

from analysis.paths import fourteenB_HI_data_path, paper1_figures_path

'''
Channel plots of the rotation subtracted HI cube.

Borrowing code from @keflavich: https://github.com/keflavich/paper_w51_evla/blob/master/plot_codes/h77a_layers.py
'''

cube_file = fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits")

cube = SpectralCube.read(cube_file)

# Begin channel map code here
Nrows = 4
Ncols = 3

# integrate over velocities to make channel maps of a set width
vstart = 308  # channels
vend = 788
vstep = 10
split = Nrows * Ncols
# if int(vend - vstart) / vend > Nrows * Ncols:
#     warnings.warn("More slices than rows*columns!")
# if int(vend - vstart) / vend < Nrows * Ncols:
#     warnings.warn("Fewer slices than rows*columns")

all_slabs = np.arange(vstart, vend + vstep, vstep, dtype=int)

# Define the average beam
beam = average_beams(cube.beams)

# Split the 48 slabs into 4 12-panel figures
for n in range(4):
    p.figure(1, figsize=(12, 20)).clf()

    vchans = all_slabs[n * split: (n + 1) * split + 1]

    fig, ax = p.subplots(Nrows, Ncols,
                         sharex=True,
                         sharey=True, num=1)

    # Convert to K km/s, but beam equivalency doesn't allow it. Doing it ad hoc
    layers = \
        [cube[start:end].moment0().value *
         beam.jtok(1.414 * u.GHz) / 1000. * u.km / u.s
         for start, end in zip(vchans[:-1], vchans[1:])]
    # Determine the maximum value to display
    mx = np.max([np.nanmax(x).value for x in layers])

    for ii, (v1, v2) in enumerate(zip(vchans[:-1], vchans[1:])):
        layer = layers[ii]

        y, x = np.unravel_index(ii, (Nrows, Ncols))

        im = ax[y, x].imshow(layer.value, origin='lower',
                             norm=ImageNormalize(vmin=-0.001,
                                                 vmax=mx,
                                                 stretch=AsinhStretch()),
                             cmap=p.cm.gray_r)
        vel1 = cube.spectral_axis.to(u.km / u.s)[v1].value
        vel2 = cube.spectral_axis.to(u.km / u.s)[v2].value
        ax[y, x].annotate("${0:.0f}<v<{1:.0f}$".format(vel2, vel1),
                          (0.53, 0.9),
                          xycoords='axes fraction', color='k',
                          fontsize=15.5)

    # fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([1.00, 0.05, 0.04, 0.9])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("K km s$^{-1}$")

    p.subplots_adjust(hspace=0,
                      wspace=0)

    for i in range(Nrows):
        for j in range(Ncols):
            if i == 0:
                ax[i, j].xaxis.set_ticks_position('top')
                p.setp(ax[i, j].get_xticklabels(), visible=False)
                ax[i, j].xaxis.set_ticklabels([])
            elif i == Nrows - 1:
                ax[i, j].xaxis.set_ticks_position('bottom')
                p.setp(ax[i, j].get_xticklabels(), visible=True)
            else:
                ax[i, j].xaxis.set_ticks_position('none')
                p.setp(ax[i, j].get_xticklabels(), visible=False)
                ax[i, j].xaxis.set_ticklabels([])

            if j == 0:
                ax[i, j].yaxis.set_ticks_position('left')
            elif j == Ncols - 1:
                ax[i, j].yaxis.set_ticks_position('right')
                p.setp(ax[i, j].get_yticklabels(), visible=False)
                ax[i, j].yaxis.set_ticklabels([])
            else:
                ax[i, j].yaxis.set_ticks_position('none')
                p.setp(ax[i, j].get_yticklabels(), visible=False)
                ax[i, j].yaxis.set_ticklabels([])

    p.subplots_adjust(hspace=0,
                      wspace=0)

    # p.draw()

    # raw_input("Next plot?")
    p.savefig(paper1_figures_path("M33_rotsub_channels_{}.pdf".format(n)),
              bbox_inches="tight", dpi=300)
    p.savefig(paper1_figures_path("M33_rotsub_channels_{}.png".format(n)),
              bbox_inches="tight", dpi=300)
    p.clf()
