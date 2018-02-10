

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
import scipy.ndimage as nd
from astropy.modeling import models, fitting
from astropy.table import Table

# from cube_analysis.masking import ppv_connectivity_masking

from paths import (fourteenB_HI_data_wGBT_path,
                   allfigs_path, alltables_path)
from constants import hi_freq, hi_mass_conversion_Jy, distance
from plotting_styles import (default_figure, twocolumn_twopanel_figure,
                             onecolumn_twopanel_figure)


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

plt.close()

# Now focus on individual features

# The noise level in one pixel in the original cube is 2.8 K. These are
# averaged over ~5 channels, so divide that by sqrt(5). In the moment maps,
# the error is the sum over the number of channels in each slab. Then multiply
# by the number of pixels in each mask.
sigma_noise = 2.8 / np.sqrt(5) * u.K
sigma_noise_Jybm = sigma_noise.to(u.Jy, beam.jtok_equiv(hi_freq))

# Linear feature in south of galaxy
rotsub_slab1 = rotsub_cube.spectral_slab(-150 * u.km / u.s, -60 * u.km / u.s)
rotsub_slab1 = rotsub_slab1.with_mask(rotsub_slab1 > sigma_noise_Jybm)

# Make a signal mask
# masked_slab1, mask_slab1 = \
#     ppv_connectivity_masking(rotsub_slab1, smooth_chans=5, min_chan=5,
#                              peak_snr=3.,
#                              min_snr=1.5, edge_thresh=1, show_plots=False,
#                              noise_map=None, verbose=True, spatial_kernel=None)

# mom0_linfeat = masked_slab1.moment0()
# Convert to K km /s from Jy/bm m / s
conv_factor = rotsub_slab.beam.jtok(hi_freq) / 1000. * u.km / u.s

mom0_linfeat = rotsub_slab1.moment0() * conv_factor.value
mom0_linfeat._unit = u.K * u.km / u.s

# Northern blob
rotsub_slab2 = rotsub_cube.spectral_slab(52 * u.km / u.s, 140 * u.km / u.s)
rotsub_slab2 = rotsub_slab2.with_mask(rotsub_slab2 > sigma_noise_Jybm)

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
fig1.show_grayscale(stretch='linear', invert=True, vmin=0, vmax=140)
fig1.show_colorbar()
fig1.ticks.set_yspacing(10 / 60.)
fig1.ticks.set_xspacing(15 / 60.)
fig1.tick_labels.set_xformat('hh:mm')
fig1.tick_labels.set_yformat('dd:mm')
fig1.hide_axis_labels()
fig1.colorbar.set_ticks([0, 20, 60, 100, 140])
fig1.colorbar.set_font(size=11)
fig1.add_label(0.05, 0.05, '-150<v<-60 km/s ', relative=True,
               bbox={"boxstyle": "square", "facecolor": "w"},
               size=10, horizontalalignment='left',
               verticalalignment='bottom')

fig2 = FITSFigure(mom0_ncloud[750:1500, 100:1000].hdu, figure=fig,
                  subplot=(1, 2, 2))
fig2.show_grayscale(stretch='linear', invert=True, vmin=0, vmax=140)
fig2.show_colorbar()
fig2.axis_labels.hide_y()
fig2.ticks.set_yspacing(10 / 60.)
fig2.ticks.set_xspacing(15 / 60.)
fig2.tick_labels.set_xformat('hh:mm')
fig2.tick_labels.set_yformat('dd:mm')
fig2.hide_axis_labels()
fig2.colorbar.set_axis_label_text('Integrated Intensity (K km/s)')
fig2.colorbar.set_axis_label_font(size=11, weight='normal')
fig2.colorbar.set_ticks([0, 20, 60, 100, 140])
fig2.colorbar.set_font(size=11)
fig2.add_label(0.05, 0.05, '50<v<140 km/s ', relative=True,
               bbox={"boxstyle": "square", "facecolor": "w"},
               size=10, horizontalalignment='left',
               verticalalignment='bottom')

plt.tight_layout()

fig.savefig(allfigs_path("HI_maps/M33_offrot_mom0_maps.pdf"),
            bbox_inches="tight")
fig.savefig(allfigs_path("HI_maps/M33_offrot_mom0_maps.png"),
            bbox_inches="tight")
plt.close()

# default_figure()
chan_width = np.abs(np.diff(rotsub_slab1.spectral_axis[:2])[0]).to(u.km / u.s)

linfeat_spatmask = rotsub_slab1.mask.include().sum(0)[100:850, 300:1200]

linfeat_mask = mom0_linfeat[100:850, 300:1200].value > \
    (linfeat_spatmask * sigma_noise_Jybm * chan_width *
     beam.jtok(hi_freq)).value
linfeat_mask = np.logical_and(linfeat_spatmask > 30, linfeat_mask)

# Require things in the mask be at least the beam size
beam_kernel = beam.as_tophat_kernel(mom0_ncloud.header['CDELT2']).array > 0

linfeat_mask = nd.binary_opening(linfeat_mask, beam_kernel)
linfeat_mask = nd.binary_closing(linfeat_mask, beam_kernel)

pix_area_sr = ((rotsub_cube.header['CDELT2'] * u.deg)**2).to(u.sr)
beams_per_pix = pix_area_sr / beam.sr

linfeat_intens = np.nansum(mom0_linfeat[100:850, 300:1200][linfeat_mask]) / \
    beam.jtok(hi_freq) * u.Jy

linfeat_mass = linfeat_intens * hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

linfeat_mass_stddev = (sigma_noise_Jybm * 1. * u.km / u.s * linfeat_mask *
                       linfeat_spatmask).sum() * \
    hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

ncloud_spatmask = rotsub_slab2.mask.include().sum(0)[750:1500, 100:1000]

ncloud_mask = mom0_ncloud[750:1500, 100:1000].value > \
    (ncloud_spatmask * sigma_noise_Jybm * chan_width *
     beam.jtok(hi_freq)).value
ncloud_mask = np.logical_and(ncloud_spatmask > 30, ncloud_mask)


# Require things in the mask be at least the beam size
ncloud_mask = nd.binary_opening(ncloud_mask, beam_kernel)
ncloud_mask = nd.binary_closing(ncloud_mask, beam_kernel)

ncloud_intens = \
    np.nansum(mom0_ncloud[750:1500, 100:1000][ncloud_mask])

ncloud_intens = np.nansum(mom0_ncloud[750:1500, 100:1000][ncloud_mask]) / \
    beam.jtok(hi_freq) * u.Jy

ncloud_mass = ncloud_intens * hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

ncloud_mass_stddev = (sigma_noise_Jybm * 1. * u.km / u.s * ncloud_mask *
                      ncloud_spatmask).sum() * \
    hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix


print("Blue-shifted mass: {0}+/-{1}"
      .format(linfeat_mass.value, linfeat_mass_stddev))
print("Red-shifted mass: {0}+/-{1}"
      .format(ncloud_mass.value, ncloud_mass_stddev))
# Blue-shifted mass: 4983498.22324+/-2132314.62835 solMass
# Red-shifted mass: 7819298.53279+/-3304023.51807 solMass


# Examine the N plume feature in more depth

# Load in the region file
import pyregion
nplume_reg = pyregion.open(fourteenB_HI_data_wGBT_path("downsamp_1kms/nplume.reg"))

nplume_slab = rotsub_slab.subcube_from_ds9region(nplume_reg)

mean_spec = nplume_slab.mean(axis=(1, 2)) * nplume_slab.beam.jtok(hi_freq) / u.Jy
total_spec = nplume_slab.sum(axis=(1, 2))

# How much mass is in this cloud? Try fitting with two Gaussians

g_HI_init = models.Gaussian1D(amplitude=40., mean=-40000,
                              stddev=10000) +  \
    models.Gaussian1D(amplitude=65, mean=0,
                      stddev=8000)

fit_g = fitting.LevMarLSQFitter()

vels = rotsub_slab.spectral_axis

g_HI = fit_g(g_HI_init, vels, total_spec)

# The covariance matrix is hidden away... tricksy
cov = fit_g.fit_info['param_cov']
parnames = g_HI.param_names
parvals = g_HI.parameters
parerrs = np.sqrt(np.diag(cov))

# In Jy/bm km/s
nplume_intens = np.sqrt(2 * np.pi) * parvals[0] * (parvals[2] / 1000.)
nplume_err = np.sqrt(2 * np.pi) * \
    (parvals[2] * parerrs[0] + parvals[0] * parerrs[2]) / 1000.

nplume_mass = nplume_intens * hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

nplume_mass_stddev = nplume_err * hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

print("N Plume mass {0}+/-{1}".format(nplume_mass.value, nplume_mass_stddev))
# N Plume mass 6455308.85877+/-312563.675469 solMass

# Now fit the mean spectrum for the figure

g_HI_init = models.Gaussian1D(amplitude=10., mean=-40000,
                              stddev=10000) +  \
    models.Gaussian1D(amplitude=16., mean=0,
                      stddev=8000)

fit_g = fitting.LevMarLSQFitter()

vels = rotsub_slab.spectral_axis

g_HI = fit_g(g_HI_init, vels, mean_spec)

cov = fit_g.fit_info['param_cov']
parnames = g_HI.param_names
parvals = g_HI.parameters
parerrs = np.sqrt(np.diag(cov))

# Make a nice figure
onecolumn_twopanel_figure()

rotsub_slab3 = rotsub_slab.spectral_slab(-64000 * u.m / u.s,
                                         -21000 * u.m / u.s)
rotsub_slab3 = rotsub_slab3[:, 1250:1550, 450:750]

mom0_nplume = rotsub_slab3.moment0() * conv_factor.value
mom0_nplume._unit = u.K * u.km / u.s

fig = plt.figure()

fig1 = FITSFigure(mom0_nplume.hdu, figure=fig,
                  subplot=(2, 1, 1))
fig1.show_grayscale(stretch='linear', invert=True)
fig1.show_colorbar()
fig1.ticks.set_xspacing(6 / 60.)
fig1.tick_labels.set_xformat('hh:mm:ss')
fig1.tick_labels.set_yformat('dd:mm:ss')
fig1.colorbar.set_axis_label_text('Integrated Intensity (K km/s)')
fig1.colorbar.set_font(size=11)
fig1.show_regions(nplume_reg)
fig1.add_label(0.05, 0.05, '-64<v<-21 km/s ', relative=True,
               bbox={"boxstyle": "square", "facecolor": "w"},
               size=10, horizontalalignment='left',
               verticalalignment='bottom')

ax = fig.add_subplot(212)

ax.plot(vels / 1000., mean_spec.value, drawstyle='steps-mid')
ax.plot(vels / 1000., g_HI[0](vels), '--', label='Cloud Comp.')
ax.plot(vels / 1000., g_HI[1](vels), '-.', label='Disk Comp.')
ax.plot(vels / 1000., g_HI(vels), '-', alpha=0.8, label='Total Fit')
ax.set_xlim([-100, 100])

ax.set_xlabel("Velocity (km/s)")
ax.set_ylabel("Intensity (K)")

ax.legend(frameon=True, loc='upper left')
ax.grid()

plt.tight_layout()

fig.savefig(allfigs_path("HI_maps/M33_nplume_mom0_wfit.pdf"),
            bbox_inches="tight")
fig.savefig(allfigs_path("HI_maps/M33_nplume_mom0_wfit.png"),
            bbox_inches="tight")
plt.close()

# Save the fit parameters and the masses

fit_tab = Table({"Names": parnames, "Values": parvals, "Errors": parerrs})

fit_tab.write(alltables_path("nplume_gaussian_fit.csv"))
