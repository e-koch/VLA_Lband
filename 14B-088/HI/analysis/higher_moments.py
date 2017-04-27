
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.cube_utils import average_beams
import numpy as np
import matplotlib.pyplot as p
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import os
from reproject import reproject_interp
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from corner import hist2d

from analysis.paths import (fourteenB_HI_data_path, paper1_figures_path,
                            iram_co21_data_path, data_path)
from constants import (hi_freq, cube_name, moment0_name, lwidth_name,
                       skew_name, kurt_name, mask_name, moment1_name,
                       peaktemps_name, peakvels_name,
                       rotsub_cube_name, rotsub_mask_name)
from plotting_styles import (twocolumn_figure, onecolumn_figure,
                             default_figure)

'''
Investigating skewness and kurtosis in the 14B-088 cube.
'''

cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))
mask = fits.open(fourteenB_HI_data_path(mask_name))[0].data > 0
cube = cube.with_mask(mask)

# Checked whether these variations in the spectra are driven by rotation
# There is not significant change on these scales.
# cube = SpectralCube.read(fourteenB_HI_data_path(rotsub_cube_name))
# mask = fits.open(fourteenB_HI_data_path(rotsub_mask_name))[0].data > 0
# cube = cube.with_mask(mask)

mom0_hdu = fits.open(fourteenB_HI_data_path(moment0_name))[0]
mom0 = Projection.from_hdu(mom0_hdu)

mom1_hdu = fits.open(fourteenB_HI_data_path(moment1_name))[0]
mom1 = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_HI_data_path(lwidth_name))[0]
lwidth = Projection.from_hdu(lwidth_hdu)

skew_hdu = fits.open(fourteenB_HI_data_path(skew_name))[0]
skew = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_HI_data_path(kurt_name))[0]
kurt = Projection.from_hdu(kurt_hdu)

peaktemps_hdu = fits.open(fourteenB_HI_data_path(peaktemps_name))[0]
peaktemps = Projection.from_hdu(peaktemps_hdu)

peakvels_hdu = fits.open(fourteenB_HI_data_path(peakvels_name))[0]
peakvels = Projection.from_hdu(peakvels_hdu)

# CO and 3.6 um for comparison
co_hdu = fits.open(iram_co21_data_path("m33.ico.hireprojection.fits"))[0]
co_map = Projection(co_hdu.data, wcs=WCS(co_hdu.header), unit=u.K)

irac1_hdu = fits.open(os.path.join(data_path, "Spitzer/irac1_3.6um/ch1_122_bgsub.fits"))[0]
irac1_reproj = reproject_interp(irac1_hdu, mom0.header)[0]
irac1_map = Projection(irac1_reproj, wcs=mom0.wcs, unit=u.MJy / u.sr)

# Make a nice 2 panel figure
twocolumn_figure(fig_ratio=0.5)
ax = p.subplot(121, projection=skew.wcs)
im = ax.imshow(skew.value,
               origin='lower', vmin=-2, vmax=2,
               interpolation='nearest', cmap='seismic')
# ax.set_title("Skewness")
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
cbar = p.colorbar(im)
cbar.set_label("Skewness")

ax2 = p.subplot(122, projection=kurt.wcs)
im2 = ax2.imshow(kurt.value, vmin=-3, vmax=3,
                 origin='lower', interpolation='nearest', cmap='seismic')
ax2.set_xlabel("RA (J2000)")
lat = ax2.coords[1]
lat.set_ticklabel_visible(False)
cbar2 = p.colorbar(im2)
cbar2.set_label("Kurtosis")

p.tight_layout()

p.savefig(paper1_figures_path("skew_kurt_maps.pdf"))
p.savefig(paper1_figures_path("skew_kurt_maps.png"))

raw_input("Next plot?")
p.close()


# Interesting regions:

# Northern HI infall: [1373, 630], [1373, 640], [1373, 650]
nplume_positions = \
    np.array([[1375, 603], [1375, 613], [1375, 623], [1375, 633], [1375, 643]])

# SArm positions
# sarm_positions = np.array([[690, 645], [684, 645], [678, 645], [672, 645],
#                            [666, 645], [660, 645], [654, 645]])
sarm_positions = np.array([[680, 720], [670, 720], [660, 720], [650, 720],
                           [640, 720],
                           [630, 720]])

# NGC 604
ngc604_positions = np.array([[980, 445], [970, 445], [960, 445], [950, 445],
                             [940, 445], [930, 445], [920, 445]])


# Southern Arm, NGC 604, Northern plume
slicer = [(slice(594, 720), slice(670, 790)),
          (slice(878, 1000), slice(395, 506)),
          (slice(1250, 1470), slice(500, 740))]

spec_posns = [sarm_positions, ngc604_positions, nplume_positions]

names = ["sarm", "ngc604", "nplume"]

twocolumn_figure(fig_ratio=0.67)

for i, (slices, posns) in enumerate(zip(slicer, spec_posns)):
    # Moment 0 of the last few channels to highlight in-falling HI plume
    if i == 2:
        loc_mom0 = cube[950:].moment0()
    else:
        loc_mom0 = mom0
    # Moment 0
    ax = p.subplot(231, projection=mom0[slices].wcs)
    ax.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
            color="chartreuse",
            marker="D", linestyle="None",
            markersize=6)
    ax.imshow(loc_mom0[slices].value, origin='lower',
              interpolation='nearest', cmap='binary')
    ax.set_ylabel("DEC (J2000)")
    ax.text(10, 10, "Moment 0", bbox={"boxstyle": "square", "facecolor": "w"},
            horizontalalignment='left', verticalalignment='bottom')
    lon = ax.coords[0]
    lon.set_ticklabel_visible(False)

    # Lwidth
    ax2 = p.subplot(232, projection=lwidth[slices].wcs)
    ax2.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
             color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax2.imshow(lwidth.value[slices], origin='lower',
               interpolation='nearest', cmap='binary')
    lon = ax2.coords[0]
    lon.set_ticklabel_visible(False)
    lat = ax2.coords[1]
    lat.set_ticklabel_visible(False)
    ax2.text(10, 10, "Line Width", bbox={"boxstyle": "square", "facecolor": "w"},
             horizontalalignment='left', verticalalignment='bottom')

    # Skewness
    ax3 = p.subplot(234, projection=skew[slices].wcs)
    ax3.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
             color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax3.imshow(skew.value[slices], vmin=-2, vmax=2,
               origin='lower',
               interpolation='nearest', cmap='seismic')
    ax3.set_ylabel("DEC (J2000)")
    ax3.set_xlabel("RA (J2000)")
    ax3.text(10, 10, "Skewness", bbox={"boxstyle": "square", "facecolor": "w"},
             horizontalalignment='left', verticalalignment='bottom')

    # Kurtosis
    ax4 = p.subplot(235, projection=kurt[slices].wcs)
    ax4.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
             color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax4.imshow(kurt.value[slices], vmin=1, vmax=4,
               origin='lower', interpolation='nearest', cmap='binary')
    lat = ax4.coords[1]
    lat.set_ticklabel_visible(False)
    ax4.text(10, 10, "Kurtosis", bbox={"boxstyle": "square", "facecolor": "w"},
             horizontalalignment='left', verticalalignment='bottom')
    ax4.set_xlabel("RA (J2000)")
    # CO
    ax5 = p.subplot(233, projection=co_map[slices].wcs)
    ax5.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
             color="chartreuse",
             marker="D", linestyle="None",
             markersize=6, markeredgecolor='k')
    ax5.imshow(co_map.value[slices],
               origin='lower', interpolation='nearest', cmap='binary')
    lon = ax5.coords[0]
    lon.set_ticklabel_visible(False)
    lat = ax5.coords[1]
    lat.set_ticklabel_visible(False)
    ax5.text(10, 10, "CO(2-1) Peak Int.", bbox={"boxstyle": "square", "facecolor": "w"},
             horizontalalignment='left', verticalalignment='bottom')
    # IRAC 3.6 um
    ax6 = p.subplot(236, projection=irac1_map[slices].wcs)
    ax6.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
             color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax6.imshow(irac1_map.value[slices],
               origin='lower', interpolation='nearest', cmap='binary',
               norm=ImageNormalize(vmin=0.0,
                                   vmax=1.,
                                   stretch=AsinhStretch()))
    ax6.set_xlabel("RA (J2000)")
    lat = ax6.coords[1]
    lat.set_ticklabel_visible(False)
    ax6.text(10, 10, r"IRAC 3.6 $\mu$m", bbox={"boxstyle": "square", "facecolor": "w"},
             horizontalalignment='left', verticalalignment='bottom')

    p.subplots_adjust(hspace=0.02, wspace=0.02)
    p.draw()

    p.savefig(paper1_figures_path("{}_moments.png").format(names[i]))
    p.savefig(paper1_figures_path("{}_moments.pdf").format(names[i]))

    # raw_input("Next plot?")
    p.clf()

p.close()

# Zoom in on the interesting portion of the spectrum
spectral_cuts = np.array([[-72, -180], [-200, -280], [-220, -315]]) * \
    u.km / u.s

# If using with rotation-subtracted cube. This does NOT change the small-scale
# variation in the components.
# spectral_cuts = np.array([[-60, 60], [-60, 60], [-60, 60]]) * \
#     u.km / u.s


# In text, put the physical distance away from the "center" (usually the
# brightest). Each pixel in the map is 12 pc (1 pix = 3'', 1'' = 4 pc)
sarm_dists = np.array([3, 2, 1, 0, -1, -2]) * 10 * 12 * u.pc
ngc604_dists = np.arange(3, -4, -1) * 10 * 12 * u.pc
nplume_dists = np.arange(2, -3, -1) * 10 * 12 * u.pc

dists = [sarm_dists, ngc604_dists, nplume_dists]

# Now save spectra going through the interesting regions:
onecolumn_figure(fig_ratio=2.27)

for posns, cuts, name, dist in zip(spec_posns, spectral_cuts, names, dists):

    num_posns = len(posns)
    fig, axes = p.subplots(num_posns, 1, sharey=True, sharex=False)
    for i, (posn, ax) in enumerate(zip(posns, axes)):
        spec = cube.spectral_slab(cuts[0], cuts[1])[:, posn[0], posn[1]]
        spec = spec.to(u.K, equivalencies=average_beams(cube.beams).jtok_equiv(hi_freq))
        velocities = cube.spectral_slab(cuts[0], cuts[1]).spectral_axis.to(u.km / u.s).value
        ax.plot(velocities, spec.value, 'b-', drawstyle='steps-mid')
        # ax.axvline(mom1[posn[0], posn[1]].to(u.km / u.s).value)
        if i < num_posns - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Velocity (km/s)")
            for label in ax.get_xticklabels()[1::2]:
                label.set_visible(False)

    for i, (posn, ax, d) in enumerate(zip(posns, axes, dist)):
        ax.text(cuts[0].value - (cuts[0] - cuts[1]).value * 0.2,
                ax.get_ybound()[1] * 0.7,
                "{0} {1}".format(int(d.value), d.unit.to_string()),
                verticalalignment='bottom')

        ax.text(cuts[0].value - (cuts[0] - cuts[1]).value * 0.95,
                ax.get_ybound()[1] * 0.55,
                "Skew={0}\n Kurt={1}".format(round(skew.value[posn[0], posn[1]], 3),
                                             round(kurt.value[posn[0], posn[1]] - 3, 3)),
                verticalalignment='bottom')

        # ax.text(cuts[0].value - (cuts[0] - cuts[1]).value * 0.95,
        #         ax.get_ybound()[1] * 0.45,
        #         "Skew={0}\n Kurt={1}\nLW={2}\nMean={3}"
        #         .format(round(skew.value[posn[0], posn[1]], 3),
        #                 round(kurt.value[posn[0], posn[1]] - 3, 3),
        #                 round(lwidth.value[posn[0], posn[1]]/1000., 3),
        #                 round(mom1.value[posn[0], posn[1]]/1000., 3)),
        #         verticalalignment='bottom')

        if i == 0:
            for label in ax.get_yticklabels()[1::2]:
                label.set_visible(False)
        else:
            for label in ax.get_yticklabels():
                label.set_visible(False)

    fig.text(0.04, 0.5, 'Intensity (K)', va='center', rotation='vertical')
    p.draw()

    # raw_input("Continue ?")

    p.savefig(paper1_figures_path("{}_spectra.png").format(name))
    p.savefig(paper1_figures_path("{}_spectra.pdf").format(name))

    p.close()


# Compare the population of values to different properties
# Convert zeroth moment to K km/s from Jy m/s
mom0_vals = (mom0.array.flatten() * mom0.beam.jtok(hi_freq) *
             u.km / u.s / 1000.).value
mom1_vals = mom1.to(u.km / u.s).array.flatten()
lwidth_vals = lwidth.to(u.km / u.s).array.flatten()
skew_vals = skew.array.flatten()
kurt_vals = kurt.array.flatten()
peaktemps_vals = peaktemps.array.flatten()
peakvels_vals = peakvels.to(u.km / u.s).array.flatten()

good_vals = np.isfinite(skew_vals)

mask_summed = cube.mask.include().sum(0)

onecolumn_figure(font_scale=1.1)
# Skew vs. kurt
hist2d(skew_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-3, 3), (-3, 6)])
p.xlabel("Skewness")
p.ylabel("Kurtosis")
p.tight_layout()

p.savefig(paper1_figures_path("skewness_vs_kurtosis.pdf"))
p.savefig(paper1_figures_path("skewness_vs_kurtosis.png"))
p.close()

# Versus integrated intensity
hist2d(mom0_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_vals)), (-3, 3)])
p.xlabel(r"Integrated Intensity (K km s$^{-1}$)")
p.ylabel("Skewness")
p.tight_layout()

p.savefig(paper1_figures_path("mom0_vs_skewness.pdf"))
p.savefig(paper1_figures_path("mom0_vs_skewness.png"))
p.close()

hist2d(mom0_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_vals)), (-3, 6)])
p.xlabel(r"Integrated Intensity (K km s$^{-1}$)")
p.ylabel("Kurtosis")
p.tight_layout()

p.savefig(paper1_figures_path("mom0_vs_kurtosis.pdf"))
p.savefig(paper1_figures_path("mom0_vs_kurtosis.png"))
p.close()

# Take a bin of summed mask values and see where this lies on the mom0 vs
# kurtosis plane
hist2d(mom0_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_vals)), (-3, 6)])
p.xlabel(r"Integrated Intensity (K km s$^{-1}$)")
p.ylabel("Kurtosis")

mask_summed_300 = mask_summed.flatten()[good_vals] == 300
p.plot(mom0_vals[good_vals][mask_summed_300],
       kurt_vals[good_vals][mask_summed_300], 'D')
p.tight_layout()
p.savefig(paper1_figures_path("mom0_vs_kurtosis_w_maskwidth300.pdf"))
p.savefig(paper1_figures_path("mom0_vs_kurtosis_w_maskwidth300.png"))
p.close()

# Compare velocity surfaces to skew and kurt.
hist2d(mom1_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mom1_vals), np.nanmax(mom1_vals)), (-3, 3)])
p.xlabel(r"Centroid Velocity (km s$^{-1}$)")
p.ylabel("Skewness")
p.tight_layout()

p.savefig(paper1_figures_path("mom1_vs_skewness.pdf"))
p.savefig(paper1_figures_path("mom1_vs_skewness.png"))
p.close()

hist2d(mom1_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mom1_vals), np.nanmax(mom1_vals)), (-3, 6)])
p.xlabel(r"Centroid Velocity (km s$^{-1}$)")
p.ylabel("Kurtosis")
p.tight_layout()

p.savefig(paper1_figures_path("mom1_vs_kurtosis.pdf"))
p.savefig(paper1_figures_path("mom1_vs_kurtosis.png"))
p.close()

hist2d(peakvels_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peakvels_vals), np.nanmax(peakvels_vals)), (-3, 3)])
p.xlabel(r"Peak Velocity (km s$^{-1}$)")
p.ylabel("Skewness")
p.tight_layout()

p.savefig(paper1_figures_path("peakvel_vs_skewness.pdf"))
p.savefig(paper1_figures_path("peakvel_vs_skewness.png"))
p.close()

hist2d(peakvels_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peakvels_vals), np.nanmax(peakvels_vals)), (-3, 6)])
p.xlabel(r"Peak Velocity (km s$^{-1}$)")
p.ylabel("Kurtosis")
p.tight_layout()

p.savefig(paper1_figures_path("peakvel_vs_kurtosis.pdf"))
p.savefig(paper1_figures_path("peakvel_vs_kurtosis.png"))
p.close()

# Versus peak temperature
hist2d(peaktemps_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps).value, np.nanmax(peaktemps).value),
              (-3, 3)])
p.xlabel("Peak Temperature (K)")
p.ylabel("Skewness")
p.tight_layout()

p.savefig(paper1_figures_path("peaktemp_vs_skewness.pdf"))
p.savefig(paper1_figures_path("peaktemp_vs_skewness.png"))
p.close()

hist2d(peaktemps_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps).value, np.nanmax(peaktemps).value),
              (-3, 6)])
p.xlabel("Peak Temperature (K)")
p.ylabel("Kurtosis")
p.tight_layout()

p.savefig(paper1_figures_path("peaktemp_vs_kurtosis.pdf"))
p.savefig(paper1_figures_path("peaktemp_vs_kurtosis.png"))
p.close()

hist2d(peaktemps_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps).value, np.nanmax(peaktemps).value),
              (-3, 6)])
p.xlabel("Peak Temperature (K)")
p.ylabel("Kurtosis")

p.plot(peaktemps_vals[good_vals][mask_summed_300],
       kurt_vals[good_vals][mask_summed_300], 'D')
p.tight_layout()
p.savefig(paper1_figures_path("peaktemp_vs_kurtosis_w_maskwidth300.pdf"))
p.savefig(paper1_figures_path("peaktemp_vs_kurtosis_w_maskwidth300.png"))
p.close()


# How badly is each affected by the shape of the mask?
hist2d(mask_summed.flatten()[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)), (-3, 3)])
p.xlabel("Summed Mask (spectral pixels)")
p.ylabel("Skewness")
p.tight_layout()

p.savefig(paper1_figures_path("maskshape_vs_skewness.pdf"))
p.savefig(paper1_figures_path("maskshape_vs_skewness.png"))
p.close()

hist2d(mask_summed.flatten()[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)), (-3, 6)])
p.xlabel("Summed Mask (spectral pixels)")
p.ylabel("Kurtosis")
p.tight_layout()

p.savefig(paper1_figures_path("maskshape_vs_kurtosis.pdf"))
p.savefig(paper1_figures_path("maskshape_vs_kurtosis.png"))
p.close()

# hist2d(mask_summed.flatten()[good_vals], peaktemps_vals[good_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)),
#               (np.nanmin(peaktemps_vals), np.nanmax(peaktemps_vals))])

# hist2d(mask_summed.flatten()[good_vals], mom0_vals[good_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)),
#               (np.nanmin(mom0).value, np.nanmax(mom0).value)])

# hist2d(mom0_vals[good_vals], peaktemps_vals[good_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mom0_vals), np.nanmax(mom0_vals)),
#               (np.nanmin(peaktemps_vals), np.nanmax(peaktemps_vals))])


default_figure()
