

'''
Figures and analysis highlighting anomalous gas features
'''

from spectral_cube import SpectralCube
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from aplpy import FITSFigure
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# from cube_analysis.masking import ppv_connectivity_masking

from paths import (fourteenB_HI_data_wGBT_path,
                   allfigs_path)
from constants import hi_freq
from plotting_styles import default_figure, twocolumn_twopanel_figure


cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.1kms.fits"))
rotsub_cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.rotation_corrected.1kms.fits"))


# Make 5 km/s channel maps on the rot sub cube.

rotsub_slab = rotsub_cube.spectral_slab(-150 * u.km / u.s, 150 * u.km / u.s)

del_v = 5 * u.km / u.s
vel_bins = np.arange(-150, 155, del_v.value, dtype=int)

# Begin channel map code here
Nrows = 4
Ncols = 3

split = Nrows * Ncols

# Define the average beam
beam = cube.beam

# Split the 48 slabs into 4 12-panel figures
for n in range(5):
    plt.figure(1, figsize=(8, 11)).clf()

    vchans = vel_bins[n * split: (n + 1) * split + 1]

    fig, ax = plt.subplots(Nrows, Ncols,
                           sharex=True,
                           sharey=True, num=1)

    # Convert to K km/s, but beam equivalency doesn't allow it. Doing it ad hoc
    layers = \
        [rotsub_cube.spectral_slab(start * u.km / u.s, end * u.km / u.s).moment0().value *
         beam.jtok(hi_freq) / 1000. * u.km / u.s
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
                             cmap=plt.cm.gray_r)
        # vel1 = cube.spectral_axis.to(u.km / u.s)[v1].value
        # vel2 = cube.spectral_axis.to(u.km / u.s)[v2].value
        ax[y, x].annotate("${0:.0f}<v<{1:.0f}$".format(v1, v2),
                          (0.53, 0.9),
                          xycoords='axes fraction', color='k',
                          fontsize=8)

    # fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([1.00, 0.05, 0.04, 0.9])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("K km s$^{-1}$")

    plt.subplots_adjust(hspace=0,
                        wspace=0)

    for i in range(Nrows):
        for j in range(Ncols):
            if i == 0:
                ax[i, j].xaxis.set_ticks_position('top')
                plt.setp(ax[i, j].get_xticklabels(), visible=False)
                ax[i, j].xaxis.set_ticklabels([])
            elif i == Nrows - 1:
                ax[i, j].xaxis.set_ticks_position('bottom')
                plt.setp(ax[i, j].get_xticklabels(), visible=True)
            else:
                ax[i, j].xaxis.set_ticks_position('none')
                plt.setp(ax[i, j].get_xticklabels(), visible=False)
                ax[i, j].xaxis.set_ticklabels([])

            if j == 0:
                ax[i, j].yaxis.set_ticks_position('left')
            elif j == Ncols - 1:
                ax[i, j].yaxis.set_ticks_position('right')
                plt.setp(ax[i, j].get_yticklabels(), visible=False)
                ax[i, j].yaxis.set_ticklabels([])
            else:
                ax[i, j].yaxis.set_ticks_position('none')
                plt.setp(ax[i, j].get_yticklabels(), visible=False)
                ax[i, j].yaxis.set_ticklabels([])

    plt.subplots_adjust(hspace=0,
                        wspace=0)

    plt.draw()

    # raw_input("Next plot?")
    plt.savefig(allfigs_path("HI_maps/M33_rotsub_channels_{}.pdf".format(n)),
                bbox_inches="tight")
    plt.savefig(allfigs_path("HI_maps/M33_rotsub_channels_{}.png".format(n)),
                bbox_inches="tight")
    plt.clf()

# Now focus on individual features

# Linear feature in south of galaxy
rotsub_slab1 = rotsub_cube.spectral_slab(-150 * u.km / u.s, -60 * u.km / u.s)

# Make a signal mask
# masked_slab1, mask_slab1 = \
#     ppv_connectivity_masking(rotsub_slab1, smooth_chans=5, min_chan=5,
#                              peak_snr=3.,
#                              min_snr=1.5, edge_thresh=1, show_plots=False,
#                              noise_map=None, verbose=True, spatial_kernel=None)

# mom0_linfeat = masked_slab1.moment0()
# Convert to K km /s from Jy/bm m / s
conv_factor = rotsub_slab1.beam.jtok(hi_freq) / 1000. * u.km / u.s

mom0_linfeat = rotsub_slab1.moment0() * conv_factor.value
mom0_linfeat._unit = u.K * u.km / u.s

# Northern blob
rotsub_slab2 = rotsub_cube.spectral_slab(52 * u.km / u.s, 140 * u.km / u.s)

# masked_slab2 = \
#     ppv_connectivity_masking(rotsub_slab2, smooth_chans=5, min_chan=5,
#                              peak_snr=3.,
#                              min_snr=1.5, edge_thresh=1, show_plots=False,
#                              noise_map=None, verbose=True)
# mom0_ncloud = masked_slab2.moment0()
mom0_ncloud = rotsub_slab2.moment0() * conv_factor.value
mom0_ncloud._unit = u.K * u.km / u.s

twocolumn_twopanel_figure()

fig = plt.figure()

fig1 = FITSFigure(mom0_linfeat[100:850, 300:1200].hdu, figure=fig,
                  subplot=(1, 2, 1))
fig1.show_grayscale(stretch='arcsinh', invert=True, vmin=0)
fig1.show_colorbar()
fig1.ticks.set_yspacing(10 / 60.)
fig1.ticks.set_xspacing(15 / 60.)
fig1.tick_labels.set_xformat('hh:mm:ss')
fig1.tick_labels.set_yformat('dd:mm:ss')
fig1.colorbar.set_ticks([0, 20, 60, 100, 140])
fig1.colorbar.set_font(size=11)
fig1.add_label(0.3, 0.1, '-150<v<-60 km/s ', relative=True,
               bbox={"boxstyle": "square", "facecolor": "w"},
               size=10)

fig2 = FITSFigure(mom0_ncloud[750:1500, 100:1000].hdu, figure=fig,
                  subplot=(1, 2, 2))
fig2.show_grayscale(stretch='arcsinh', invert=True, vmin=0)
fig2.show_colorbar()
fig2.axis_labels.hide_y()
fig2.ticks.set_yspacing(10 / 60.)
fig2.ticks.set_xspacing(15 / 60.)
fig2.tick_labels.set_xformat('hh:mm:ss')
fig2.tick_labels.set_yformat('dd:mm:ss')
fig2.colorbar.set_axis_label_text('Integrated Intensity (K km/s)')
fig2.colorbar.set_axis_label_font(size=11, weight='normal')
fig2.colorbar.set_ticks([0, 20, 60, 100, 140])
fig2.colorbar.set_font(size=11)
fig2.add_label(0.75, 0.9, '50<v<140 km/s ', relative=True,
               bbox={"boxstyle": "square", "facecolor": "w"},
               size=10)

plt.tight_layout()

fig.savefig(allfigs_path("HI_maps/M33_offrot_mom0_maps.pdf"),
            bbox_inches="tight")
fig.savefig(allfigs_path("HI_maps/M33_offrot_mom0_maps.png"),
            bbox_inches="tight")

default_figure()
