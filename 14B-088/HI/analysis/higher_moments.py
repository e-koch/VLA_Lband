
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.cube_utils import average_beams
import numpy as np
import matplotlib.pyplot as p
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from analysis.paths import fourteenB_HI_data_path, paper1_figures_path
from constants import (hi_freq, cube_name, moment0_name, lwidth_name,
                       skew_name, kurt_name, mask_name)

'''
Investigating skewness and kurtosis in the 14B-088 cube.
'''

cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))
mask = fits.open(fourteenB_HI_data_path(mask_name))[0].data > 0
cube = cube.with_mask(mask)

mom0_hdu = fits.open(fourteenB_HI_data_path(moment0_name))[0]
mom0 = Projection(mom0_hdu.data, wcs=WCS(mom0_hdu.header),
                  unit=u.Jy * u.m / u.s)

lwidth_hdu = fits.open(fourteenB_HI_data_path(lwidth_name))[0]
lwidth = Projection(lwidth_hdu.data, wcs=WCS(lwidth_hdu.header),
                    unit=u.m / u.s)

skew_hdu = fits.open(fourteenB_HI_data_path(skew_name))[0]
skew = Projection(skew_hdu.data, wcs=WCS(skew_hdu.header), unit=u.Unit(""))

kurt_hdu = fits.open(fourteenB_HI_data_path(kurt_name))[0]
kurt = Projection(kurt_hdu.data, wcs=WCS(kurt_hdu.header), unit=u.Unit(""))


# Make a nice 2 panel figure
ax = p.subplot(121, projection=skew.wcs)
ax.imshow(skew.value,
          origin='lower', vmin=-2, vmax=2,
          interpolation='nearest', cmap='seismic')
ax.set_title("Skewness")
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
ax2 = p.subplot(122, projection=kurt.wcs)
ax2.imshow(kurt.value, vmin=1, vmax=4,
           origin='lower', interpolation='nearest', cmap='binary')
ax2.set_title("Kurtosis")
ax2.set_xlabel("RA (J2000)")
lat = ax2.coords[1]
lat.set_ticklabel_visible(False)

p.savefig(paper1_figures_path("skew_kurt_maps.pdf"))
p.savefig(paper1_figures_path("skew_kurt_maps.png"))

# raw_input("Next plot?")
p.clf()

# Interesting regions:

# Northern HI infall: [1373, 630], [1373, 640], [1373, 650]
nplume_positions = \
    np.array([[1375, 603], [1375, 613], [1375, 623], [1375, 633], [1375, 643]])

# SArm positions
sarm_positions = np.array([[680, 730], [670, 730], [660, 730], [650, 730],
                           [640, 730],
                           [630, 730]])

# NGC 604
ngc604_positions = np.array([[980, 445], [970, 445], [960, 445], [950, 445],
                             [940, 445], [930, 445], [920, 445]])


# Southern Arm, NGC 604, Northern plume
slicer = [(slice(594, 794), slice(620, 889)),
          (slice(868, 1017), slice(365, 536)),
          (slice(1250, 1470), slice(500, 740))]

spec_posns = [sarm_positions, ngc604_positions, nplume_positions]

names = ["sarm", "ngc604", "nplume"]

for i, (slices, posns) in enumerate(zip(slicer, spec_posns)):
    # Moment 0 of the last few channels to highlight in-falling HI plume
    if i == 2:
        loc_mom0 = cube[950:].moment0()
    else:
        loc_mom0 = mom0
    # Moment 0
    ax = p.subplot(221, projection=mom0[slices].wcs)
    ax.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
            color="chartreuse",
            marker="D", linestyle="None",
            markersize=6)
    ax.imshow(loc_mom0[slices].value, origin='lower',
              interpolation='nearest', cmap='binary')
    ax.set_ylabel("DEC (J2000)")
    ax.set_title("Moment 0")
    lon = ax.coords[0]
    lon.set_ticklabel_visible(False)

    # Lwidth
    ax2 = p.subplot(222, projection=lwidth[slices].wcs)
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
    ax2.set_title("Linewidth")

    # Skewness
    ax3 = p.subplot(223, projection=skew[slices].wcs)
    ax3.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
             color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax3.imshow(skew.value[slices], vmin=-2, vmax=2,
               origin='lower',
               interpolation='nearest', cmap='seismic')
    ax3.set_ylabel("DEC (J2000)")
    ax3.set_xlabel("RA (J2000)")
    ax3.set_title("Skewness")
    # Kurtosis
    ax4 = p.subplot(224, projection=kurt[slices].wcs)
    ax4.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start,
             color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax4.imshow(kurt.value[slices], vmin=1, vmax=4,
               origin='lower', interpolation='nearest', cmap='binary')
    ax4.set_xlabel("RA (J2000)")
    lat = ax4.coords[1]
    lat.set_ticklabel_visible(False)
    ax4.set_title("Kurtosis")

    p.draw()

    p.savefig(paper1_figures_path("{}_moments.png").format(names[i]))
    p.savefig(paper1_figures_path("{}_moments.pdf").format(names[i]))

    # raw_input("Next plot?")
    p.clf()

# Zoom in on the interesting portion of the spectrum
spectral_cuts = np.array([[-72, -180], [-200, -280], [-220, -315]]) * \
    u.km / u.s

# Now save spectra going through the interesting regions:
for posns, cuts, name in zip(spec_posns, spectral_cuts, names):

    num_posns = len(posns)
    fig, axes = p.subplots(num_posns, 1, sharey=True, sharex=False,
                           figsize=(6, 10))
    for i, (posn, ax) in enumerate(zip(posns, axes)):
        spec = cube.spectral_slab(cuts[0], cuts[1])[:, posn[0], posn[1]]
        spec = spec.to(u.K, equivalencies=average_beams(cube.beams).jtok_equiv(hi_freq))
        velocities = cube.spectral_slab(cuts[0], cuts[1]).spectral_axis.to(u.km / u.s).value
        ax.plot(velocities, spec.value, 'b-', drawstyle='steps-mid')
        if i < num_posns - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Velocity (km/s)")
            for label in ax.get_xticklabels()[1::2]:
                label.set_visible(False)

    for i, ax in enumerate(axes):
        if i == 0:
            for label in ax.get_yticklabels()[1::2]:
                label.set_visible(False)
        else:
            for label in ax.get_yticklabels():
                label.set_visible(False)

    fig.text(0.04, 0.5, 'Intensity (K)', va='center', rotation='vertical')
    p.draw()

    p.savefig(paper1_figures_path("{}_spectra.png").format(name))
    p.savefig(paper1_figures_path("{}_spectra.pdf").format(name))

    # raw_input("Next plot?")
    p.close()
