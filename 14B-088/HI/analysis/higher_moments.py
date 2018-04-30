
from spectral_cube import SpectralCube, Projection
import numpy as np
import matplotlib.pyplot as p
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import os
from reproject import reproject_interp
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from paths import (fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   allfigs_path,
                   iram_co21_14B088_data_path, data_path)
from constants import hi_freq
from plotting_styles import (twocolumn_figure, onecolumn_figure,
                             default_figure)

'''
Investigating skewness and kurtosis in the 14B-088 cube.
'''

if not os.path.exists(allfigs_path("HI_maps")):
    os.mkdir(allfigs_path("HI_maps"))


cube = SpectralCube.read(fourteenB_HI_file_dict["Cube"])
mask = fits.open(fourteenB_HI_file_dict["Source_Mask"])[0].data > 0
cube = cube.with_mask(mask)

# Checked whether these variations in the spectra are driven by rotation
# There is not significant change on these scales.
# cube = SpectralCube.read(fourteenB_HI_data_path(rotsub_cube_name))
# mask = fits.open(fourteenB_HI_data_path(rotsub_mask_name))[0].data > 0
# cube = cube.with_mask(mask)

mom0_hdu = fits.open(fourteenB_HI_file_dict['Moment0'])[0]
mom0 = Projection.from_hdu(mom0_hdu)

mom1_hdu = fits.open(fourteenB_HI_file_dict['Moment1'])[0]
mom1 = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_HI_file_dict["LWidth"])[0]
lwidth = Projection.from_hdu(lwidth_hdu)

skew_hdu = fits.open(fourteenB_HI_file_dict['Skewness'])[0]
skew = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_HI_file_dict['Kurtosis'])[0]
kurt = Projection.from_hdu(kurt_hdu)

peaktemps_hdu = fits.open(fourteenB_wGBT_HI_file_dict["PeakTemp"])[0]
peaktemps = Projection.from_hdu(peaktemps_hdu)

peakvels_hdu = fits.open(fourteenB_HI_file_dict['PeakVels'])[0]
peakvels = Projection.from_hdu(peakvels_hdu)

# Feathered Skew and Kurt
skew_feather_hdu = fits.open(fourteenB_wGBT_HI_file_dict['Skewness'])[0]
skew_feather = Projection.from_hdu(skew_feather_hdu)

kurt_feather_hdu = fits.open(fourteenB_wGBT_HI_file_dict['Kurtosis'])[0]
kurt_feather = Projection.from_hdu(kurt_feather_hdu)

# CO and 3.6 um for comparison
co_hdu = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits"))[0]
co_map = Projection(co_hdu.data, wcs=WCS(co_hdu.header), unit=u.K)

irac1_hdu = fits.open(os.path.join(data_path, "Spitzer/irac1_3.6um/ch1_122_bgsub.fits"))[0]
irac1_reproj = reproject_interp(irac1_hdu, mom0.header)[0]
irac1_map = Projection(irac1_reproj, wcs=mom0.wcs, unit=u.MJy / u.sr)

# Make a nice 2 panel figure
twocolumn_figure(fig_ratio=0.5)
ax = p.subplot(121, projection=skew.wcs)
im = ax.imshow(skew.value,
               origin='lower', vmin=-3, vmax=3,
               interpolation='nearest', cmap='seismic')
# ax.set_title("Skewness")
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')

cbar = p.colorbar(im)
cbar.set_label("Skewness")

ax2 = p.subplot(122, projection=kurt.wcs)
im2 = ax2.imshow(kurt.value, vmin=-3, vmax=3,
                 origin='lower', interpolation='nearest', cmap='seismic')
ax2.set_xlabel("RA (J2000)")
lon = ax2.coords[0]
lon.set_major_formatter('hh:mm')
lat = ax2.coords[1]
lat.set_ticklabel_visible(False)
cbar2 = p.colorbar(im2)
cbar2.set_label("Kurtosis")

p.tight_layout()

p.savefig(allfigs_path("HI_maps/skew_kurt_maps.pdf"))
p.savefig(allfigs_path("HI_maps/skew_kurt_maps.png"))
p.close()

# And another with the feathered data
ax = p.subplot(121, projection=skew_feather.wcs)
im = ax.imshow(skew_feather.value,
               origin='lower', vmin=-3, vmax=3,
               interpolation='nearest', cmap='seismic')
# ax.set_title("Skewness")
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')

cbar = p.colorbar(im)
cbar.set_label("Skewness")

ax2 = p.subplot(122, projection=kurt_feather.wcs)
im2 = ax2.imshow(kurt_feather.value, vmin=-3, vmax=3,
                 origin='lower', interpolation='nearest', cmap='seismic')
ax2.set_xlabel("RA (J2000)")
lat = ax2.coords[1]
lat.set_ticklabel_visible(False)
lon = ax2.coords[0]
lon.set_major_formatter('hh:mm')

cbar2 = p.colorbar(im2)
cbar2.set_label("Kurtosis")

p.tight_layout()

p.savefig(allfigs_path("HI_maps/skew_kurt_feather_maps.pdf"))
p.savefig(allfigs_path("HI_maps/skew_kurt_feather_maps.png"))
p.close()


# Now figures of just skewness and just kurtosis

ax = p.subplot(121, projection=skew_feather.wcs)
im = ax.imshow(skew.value,
               origin='lower', vmin=-3, vmax=3,
               interpolation='nearest', cmap='seismic')
# ax.set_title("Skewness")
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')

ax2 = p.subplot(122, projection=kurt_feather.wcs)
im2 = ax2.imshow(skew_feather.value, vmin=-3, vmax=3,
                 origin='lower', interpolation='nearest', cmap='seismic')
ax2.set_xlabel("RA (J2000)")
lon = ax2.coords[0]
lon.set_major_formatter('hh:mm')

lat = ax2.coords[1]
lat.set_ticklabel_visible(False)
cbar2 = p.colorbar(im2)
cbar2.set_label("Skewness")

p.tight_layout()

p.savefig(allfigs_path("HI_maps/skew_maps.pdf"))
p.savefig(allfigs_path("HI_maps/skew_maps.png"))
p.close()

ax = p.subplot(121, projection=skew_feather.wcs)
im = ax.imshow(kurt.value,
               origin='lower', vmin=-3, vmax=3,
               interpolation='nearest', cmap='seismic')
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')

ax2 = p.subplot(122, projection=kurt_feather.wcs)
im2 = ax2.imshow(kurt_feather.value, vmin=-3, vmax=3,
                 origin='lower', interpolation='nearest', cmap='seismic')
ax2.set_xlabel("RA (J2000)")
lat = ax2.coords[1]
lat.set_ticklabel_visible(False)
lon = ax2.coords[0]
lon.set_major_formatter('hh:mm')

cbar2 = p.colorbar(im2)
cbar2.set_label("Kurtosis")

p.tight_layout()

p.savefig(allfigs_path("HI_maps/kurt_maps.pdf"))
p.savefig(allfigs_path("HI_maps/kurt_maps.png"))
p.close()

ax = p.subplot(121, projection=skew_feather.wcs)
im = ax.imshow(peaktemps.value,
               origin='lower',
               interpolation='nearest', cmap='binary')
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
lon = ax.coords[0]
lon.set_major_formatter('hh:mm')

cbar = p.colorbar(im)
cbar.set_label("Peak HI Temperature (K)")

ax2 = p.subplot(122, projection=kurt_feather.wcs)
im2 = ax2.imshow(kurt_feather.value, vmin=-3, vmax=3,
                 origin='lower', interpolation='nearest', cmap='seismic')
ax2.set_xlabel("RA (J2000)")
lat = ax2.coords[1]
lat.set_ticklabel_visible(False)
lon = ax2.coords[0]
lon.set_major_formatter('hh:mm')

cbar2 = p.colorbar(im2)
cbar2.set_label("Kurtosis")

# p.tight_layout()

p.savefig(allfigs_path("HI_maps/peaktemp_kurt_feath_maps.pdf"))
p.savefig(allfigs_path("HI_maps/peaktemp_kurt_feath_maps.png"))
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

    p.savefig(allfigs_path("HI_maps/{}_moments.png").format(names[i]))
    p.savefig(allfigs_path("HI_maps/{}_moments.pdf").format(names[i]))

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
        spec = spec.to(u.K, equivalencies=cube.beams.average_beam().jtok_equiv(hi_freq))
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

    p.savefig(allfigs_path("HI_maps/{}_spectra.png").format(name))
    p.savefig(allfigs_path("HI_maps/{}_spectra.pdf").format(name))

    p.close()

default_figure()
