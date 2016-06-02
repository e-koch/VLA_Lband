
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import numpy as np
import matplotlib.pyplot as p
import os
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

'''
Investigating skewness and kurtosis in the 14B-088 cube.
'''

save_moments = False

data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"

cube = \
    SpectralCube.read(os.path.join(data_path,
                                    "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))

mask = fits.getdata(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked_source_mask.fits"))
mask = mask.astype(bool)

cube = cube.with_mask(mask)

# Load in the linewidth
lwidth_hdu = fits.open(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.lwidth.fits"))[0]
lwidth = Projection(lwidth_hdu.data, wcs=WCS(lwidth_hdu.header), unit=u.km / u.s)

# And a moment 0
mom0_hdu = fits.open(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits"))[0]
mom0 = Projection(mom0_hdu.data, wcs=WCS(mom0_hdu.header), unit=u.Jy * u.km / u.s)


# Skewness: Converting from m^3/s^3 to km/s
if save_moments:
    mom3 = cube.moment(order=3)
    mom3.write(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom3.fits"),
               overwrite=True)
else:
    mom3_hdu = fits.open(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom3.fits"))[0]
    mom3 = Projection(mom3_hdu.data, wcs=WCS(mom3_hdu.header), unit=(u.km / u.s)**3)

# Kurtosis: Subtract 3 to center on 0 (assuming a Gaussian)
if save_moments:
    mom4 = cube.moment(order=4)
    mom4.write(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom4.fits"),
               overwrite=True)
else:
    mom4_hdu = fits.open(os.path.join(data_path, "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom4.fits"))[0]
    mom4 = Projection(mom4_hdu.data, wcs=WCS(mom4_hdu.header), unit=(u.km / u.s)**4)
# mom4_resc = (mom4.value - 3) ** 0.25 / 1000.

# Make a nice 2 panel figure
ax = p.subplot(121, projection=mom3.wcs)
ax.imshow(np.arctan(mom3.value / np.nanpercentile(mom3.value, 85)),
          origin='lower',
          interpolation='nearest', cmap='seismic')
ax.set_title("Skewness")
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
ax2 = p.subplot(122, projection=mom4.wcs)
ax2.imshow(np.arctan(mom4.value / np.nanpercentile(mom4.value, 85)),
           origin='lower', interpolation='nearest', cmap='binary')
ax2.set_title("Kurtosis")
ax2.set_xlabel("RA (J2000)")
lat = ax2.coords[1]
lat.set_ticklabel_visible(False)

raw_input("Next plot?")
p.clf()

# Interesting regions:

# Northern HI infall: [1373, 630], [1373, 640], [1373, 650]
nplume_positions = \
    np.array([[1385, 603], [1385, 613], [1385, 623], [1385, 633], [1385, 643]])

# SArm positions
sarm_positions = np.array([[670, 740], [660, 740], [650, 740], [640, 740],
                           [630, 740]])

# NGC 604
ngc604_positions = np.array([[980, 455], [970, 455], [960, 455], [950, 455],
                             [940, 455], [930, 455], [920, 455]])


# Southern Arm, NGC 604, Northern plume
slicer = [(slice(594, 794), slice(620, 889)),
          (slice(868, 1017), slice(365, 536)),
          (slice(1250, 1470), slice(500, 740))]

spec_posns = [sarm_positions, ngc604_positions, nplume_positions]

for i, (slices, posns) in enumerate(zip(slicer, spec_posns)):
    # Moment 0 of the last few channels to highlight in-falling HI plume
    if i == 2:
        loc_mom0 = cube[950:].moment0()
    else:
        loc_mom0 = mom0
    # Moment 0
    ax = p.subplot(221, projection=mom0[slices].wcs)
    ax.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start, color="chartreuse",
            marker="D", linestyle="None",
            markersize=6)
    ax.imshow(loc_mom0.value[slices], origin='lower',
              interpolation='nearest', cmap='binary')
    ax.set_ylabel("DEC (J2000)")
    ax.set_title("Moment 0")
    lon = ax.coords[0]
    lon.set_ticklabel_visible(False)

    # Lwidth
    ax2 = p.subplot(222, projection=lwidth[slices].wcs)
    ax2.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start, color="chartreuse",
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
    ax3 = p.subplot(223, projection=mom3[slices].wcs)
    ax3.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start, color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax3.imshow(np.arctan(mom3[slices].value / np.nanpercentile(mom3.value, 85)),
               origin='lower',
               interpolation='nearest', cmap='seismic')
    ax3.set_ylabel("DEC (J2000)")
    ax3.set_xlabel("RA (J2000)")
    ax3.set_title("Skewness")
    # Kurtosis
    ax4 = p.subplot(224, projection=mom4[slices].wcs)
    ax4.plot(posns[:, 1] - slices[1].start, posns[:, 0] - slices[0].start, color="chartreuse",
             marker="D", linestyle="None",
             markersize=6)
    ax4.imshow(np.arctan(mom4[slices].value / np.nanpercentile(mom4.value, 85)),
               origin='lower', interpolation='nearest', cmap='binary')
    ax4.set_xlabel("RA (J2000)")
    lat = ax4.coords[1]
    lat.set_ticklabel_visible(False)
    ax4.set_title("Kurtosis")

    p.draw()
    raw_input("Next plot?")
    p.clf()

# Zoom in on the interesting portion of the spectrum
spectral_cuts = np.array([[-72, -180], [-200, -280], [-220, -315]]) * u.km / u.s

# Now save spectra going through the interesting regions:
for posns, cuts in zip(spec_posns, spectral_cuts):

    num_posns = len(posns)
    fig, axes = p.subplots(num_posns, 1, sharey=True, sharex=False, figsize=(6, 10))
    for i, (posn, ax) in enumerate(zip(posns, axes)):
        spec = cube.spectral_slab(cuts[0], cuts[1])[:, posn[0], posn[1]].to(u.K, equivalencies=cube.beam.jtok_equiv(1.414 * u.GHz))
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

    raw_input("Next plot?")
    p.close()
