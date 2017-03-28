
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as p
from scipy import ndimage as nd
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.utils.console import ProgressBar
from astropy.wcs.utils import proj_plane_pixel_area
from astropy.io import fits
from radio_beam import Beam
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import hist
from reproject import reproject_interp
import seaborn as sb

from rotation_curves.vrot_fit import return_smooth_model
from paths import (fourteenB_HI_data_path, iram_co21_data_path,
                   paper1_figures_path)
from galaxy_params import gal
from plotting_styles import onecolumn_figure, default_figure, twocolumn_figure
from constants import moment0_name, hi_freq, moment1_name, peakvels_name

try:
    from corner import hist2d
except ImportError:
    raise ImportError("You need to install the corner package!")

'''
Is bright CO better correlated w/ the HI velocity?

Compares Tpeak_CO vs |V_peak,CO - V_rot,HI|
'''

table_name = fourteenB_HI_data_path("diskfit_noasymm_noradial_nowarp_output/rad.out.csv")
tab = Table.read(table_name)

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
del cube._header[""]

# Galactic radii
radii = gal.radius(header=cube.header)

hi_model = return_smooth_model(tab, cube.header, gal) * u.m / u.s

# Now load in the HI centroids
moment1 = fits.open(fourteenB_HI_data_path(moment1_name))[0]

# And the peak HI velocities
peakvel = fits.open(fourteenB_HI_data_path(peakvels_name))[0]


# Reproject to the CO
moment1_reproj = reproject_interp(moment1, cube.wcs.celestial,
                                  shape_out=cube.shape[1:])[0] * u.m / u.s
peakvel_reproj = reproject_interp(peakvel, cube.wcs.celestial,
                                  shape_out=cube.shape[1:])[0] * u.m / u.s


# Mask everything at the edges of the map.
r_max = 6. * u.kpc
hi_model[radii > r_max] = np.NaN * u.m / u.s
cube = cube.with_mask(radii < r_max)

# Use avg noise reported in Druard+14
co_avg_noise = 20.33 * u.mK

# Only keep peaks above 3 sigma.
min_snr = 5.
min_pts = 2

Tpeak = np.zeros(cube.shape[1:]) * u.K
vel_Tpeak = np.zeros(cube.shape[1:]) * u.m / u.s

# Compute the first moment array for the whole cube.
vel_centroid = cube.with_mask(cube > co_avg_noise * 3).moment1()
# We'll zero out any noise dominated regions below, so don't
# worry about the masking.

good_posns = np.where(cube.mask.include().sum(0) > 0)

for y, x in ProgressBar(zip(*good_posns)):
    specind = cube[:, y, x].argmax()
    Tpeak[y, x] = cube[:, y, x][specind]

    if Tpeak[y, x] < min_snr * co_avg_noise:
        Tpeak[y, x] = np.NaN * u.K
        vel_Tpeak[y, x] = np.NaN * u.m / u.s
        vel_centroid[y, x] = np.NaN * u.m / u.s
    else:
        vel_Tpeak[y, x] = cube.spectral_axis[specind]

Tpeak[Tpeak == 0.0 * u.K] = np.NaN * u.K
vel_Tpeak[vel_Tpeak == 0.0 * u.m / u.s] = np.NaN * u.m / u.s

vel_diff = np.abs(hi_model - vel_Tpeak).to(u.km / u.s)

# Clean-up the points by requiring they be part of a region larger than the
# beam
Tpeak_mask = np.isfinite(Tpeak)
pixscale = proj_plane_pixel_scales(cube.wcs)[0]
beam_kern = cube.beam.as_tophat_kernel(pixscale)

# Open and close
Tpeak_mask = nd.binary_opening(Tpeak_mask, beam_kern)
Tpeak_mask = nd.binary_closing(Tpeak_mask, beam_kern)

vel_diff_values = vel_diff.value[Tpeak_mask & np.isfinite(Tpeak)]
Tpeak_values = Tpeak.value[Tpeak_mask & np.isfinite(Tpeak)]

onecolumn_figure(font_scale=1.1)
hist2d(vel_diff_values, Tpeak_values, bins=16, data_kwargs={"alpha": 0.6},
       range=[(0.0, 30.0),
              (0.0, np.max(Tpeak_values))])
p.hlines(co_avg_noise.to(u.K).value * min_snr, 0.0, 30.0, color='r',
         linestyle='--')
p.ylabel(r"T$_\mathrm{peak, CO}$ (K)")
p.xlabel(r"|V$_\mathrm{peak, CO}$ - V$_\mathrm{rot, HI}$| (km/s)")
p.grid()
p.tight_layout()

p.savefig(paper1_figures_path("co21_Tpeak_velocity_offset.pdf"))
p.savefig(paper1_figures_path("co21_Tpeak_velocity_offset.png"))
p.close()

# HI vs CO peak velocities
vel_diff_peaks = \
    (peakvel_reproj - vel_Tpeak).to(u.km / u.s)
vel_diff_peaks_values = \
    vel_diff_peaks.value[Tpeak_mask & np.isfinite(Tpeak)]

hist2d(vel_diff_peaks_values, Tpeak_values, bins=16, data_kwargs={"alpha": 0.6},
       range=[(0.0, 30.0),
              (0.0, np.max(Tpeak_values))])
p.hlines(co_avg_noise.to(u.K).value * min_snr, 0.0, 30.0, color='r',
         linestyle='--')
p.ylabel(r"T$_\mathrm{peak, CO}$ (K)")
p.xlabel(r"|V$_\mathrm{peak, CO}$ - V$_\mathrm{peak, HI}$| (km/s)")
p.grid()
p.tight_layout()

p.savefig(paper1_figures_path("co21_Tpeak_peakvelocity_offset.pdf"))
p.savefig(paper1_figures_path("co21_Tpeak_peakvelocity_offset.png"))
p.close()

# Now do the difference between HI and CO centroids
# vel_diff = (hi_model - vel_Tpeak).to(u.km / u.s)
# vel_diff_values = vel_diff[Tpeak_mask & np.isfinite(Tpeak)]

vel_diff_centroids = \
    (moment1_reproj - vel_centroid.quantity).to(u.km / u.s)
vel_diff_centroid_values = \
    vel_diff_centroids.value[Tpeak_mask & np.isfinite(Tpeak)]

# vel_diff_hicentroid_copeak = \
#     (moment1_reproj - vel_Tpeak).to(u.km / u.s)
# vel_diff_hicentroid_copeak_values = \
#     vel_diff_hicentroid_copeak.value[Tpeak_mask & np.isfinite(Tpeak)]
# vel_diff_hicentroid_copeak_values = \
#     vel_diff_hicentroid_copeak_values[np.isfinite(vel_diff_hicentroid_copeak_values)]

# one = hist(vel_diff_values, bins='scott', color='r', normed=True,
#            label='V$_\mathrm{peak, CO}$ - V$_\mathrm{rot, HI}', histtype='step')
# two = hist(vel_diff_centroid_values, bins='scott', color='g', normed=True,
#            label='V$_\mathrm{cent, CO}$ - V$_\mathrm{cent, HI}', histtype='step')
# three = hist(vel_diff_hicentroid_copeak_values, bins='scott', color='b', normed=True,
#              label='V$_\mathrm{peak, CO}$ - V$_\mathrm{cent, HI}', histtype='step')
# p.legend(frameon=True)
# p.grid()

hist2d(np.abs(vel_diff_centroid_values), Tpeak_values, bins=12,
       data_kwargs={"alpha": 0.6},
       range=[(0.0, 30.0),
              (0.0, np.max(Tpeak_values))])
p.hlines(co_avg_noise.to(u.K).value * min_snr, 0.0, 120.0, color='r',
         linestyle='--')
p.ylabel(r"T$_\mathrm{peak, CO}$ (K)")
p.xlabel(r"|V$_\mathrm{cent, CO}$ - V$_\mathrm{cent, HI}$| (km/s)")
p.grid()
p.tight_layout()

p.savefig(paper1_figures_path("co21_Tpeak_velocity_offset_centroids.pdf"))
p.savefig(paper1_figures_path("co21_Tpeak_velocity_offset_centroids.png"))
p.close()

# Where do the outliers occur
moment0 = fits.open(fourteenB_HI_data_path(moment0_name))[0]

moment0_reproj = reproject_interp(moment0, cube.wcs.celestial,
                                  shape_out=vel_diff.shape)[0]

beam = Beam.from_fits_header(moment0.header)

moment0_Kkm_s = beam.jtok(hi_freq).value * moment0_reproj / 1000.

pixscale = np.sqrt(proj_plane_pixel_area(cube.wcs.celestial))

twocolumn_figure(fig_ratio=0.95, font_scale=1.2)
sb.set_palette("colorblind")

ax = p.subplot(111, projection=cube.wcs.celestial)
im = ax.imshow(moment0_Kkm_s,
               origin='lower',
               interpolation='nearest',
               norm=ImageNormalize(vmin=-0.001,
                                   vmax=np.nanmax(moment0_Kkm_s),
                                   stretch=AsinhStretch()))
ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")
ax.set_ylim([260, 900])
ax.set_xlim([130, 620])
# ax.add_patch(beam.ellipse_to_plot(int(0.05 * moment0.shape[0]),
#                                   int(0.05 * moment0.shape[1]), pixscale))

# Where are the velocity outliers spatially?
# Assign colours to points based on their percentile in the distribution

percs = np.percentile(vel_diff_values[vel_diff_values < 50],
                      [85, 95, 99, 99.5])
labels = ["85%", "95%", "99%", "99.5%"]
labels = ["> {0} km/s ({1})".format(round(perc, 1), label)
          for perc, label in zip(percs, labels)]

cols = ["D", "s", "o", "^"]

nonan_vel_diff = vel_diff.copy()
nonan_vel_diff[np.isnan(nonan_vel_diff)] = 0.0 * u.km / u.s
nonan_vel_diff *= Tpeak_mask * np.isfinite(Tpeak)

for perc, col, label in zip(percs, cols, labels):
    ypos, xpos = np.where(np.abs(nonan_vel_diff.value) > perc)
    p.plot(xpos, ypos, col, label=label)

p.legend(loc='lower left', frameon=True)

cbar = p.colorbar(im)
cbar.set_label(r"Integrated Intensity (K km s$^{-1}$)")

p.savefig(paper1_figures_path("co21_HI_velocity_offset_outliers.pdf"))
p.savefig(paper1_figures_path("co21_HI_velocity_offset_outliers.png"))
p.close()

default_figure()
