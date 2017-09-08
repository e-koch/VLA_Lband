
'''
Update the CO21 intensities in the dlfit.fits file, which contains fits to
the Draine models for M33. The old version used an older and incomplete IRAM
CO(2-1) map.
'''

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import astropy.units as u
from spectral_cube import Projection
from radio_beam import Beam
from os.path import join as osjoin
import os
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib as mpl
from corner import hist2d, corner
from scipy.stats import theilslopes

from paths import iram_co21_data_path, data_path, allfigs_path
from plotting_styles import default_figure


fig_path = allfigs_path("Draine_Models")

# Make figure for plots
if not os.path.exists(fig_path):
    os.mkdir(fig_path)


default_figure()

# Increase the number of colours in the palette
sb.set_palette("colorblind", 10)


def standardize(data):
    return (data - np.nanmean(data)) / np.nanstd(data)


def line_distance(x, y, m, b):
    '''
    Calculate distances from a line of slope m with intercept b.
    '''

    # dists = np.abs(b + m * x - y) / np.sqrt(1 + m**2)

    # Don't apply the abs here for a true "distance" to get which points
    # are above and below the line
    dists = (y - b - m * x) / np.sqrt(1 + m**2)

    return dists


# Use the updated version, once finished
oldtab = Table.read(osjoin(data_path, "dlfit.fits"))
tab = Table.read(osjoin(data_path, "updated_dlfit.fits"))


# Plot the Umin, CO21 vs. radius
fig, axs = plt.subplots(4, 2, sharex=True, sharey=True)

umin_co10_corr = []
umin_co21_corr = []
umin_hi_corr = []

for i, ax in zip(range(7), axs.ravel()):
    mask = (tab['RGAL'] > 1000 * i) & (tab["RGAL"] < 1000 * (i + 1))
    mask = mask & np.isfinite(tab['CO21']) & np.isfinite(tab['CO10'])
    mask = mask & np.isfinite(tab['UMIN']) & np.isfinite(tab['HI'])
    mask = mask & (tab['CO10'] > 0.0)  # There's a bunch of CO10 set to 0
    # mask = mask & (tab['CO21'] > 0.05)
    # mask = mask & (tab['UMIN'] < 2)

    stand_umin = standardize(np.array(tab['UMIN'][mask]))
    stand_co21 = standardize(np.array(tab['CO21'][mask]))
    stand_co10 = standardize(np.array(tab['CO10'][mask]))
    stand_hi = standardize(np.array(tab['HI'][mask]))

    co10_fit = theilslopes(stand_co10, x=stand_umin, alpha=0.997)
    co21_fit = theilslopes(stand_co21, x=stand_umin, alpha=0.997)
    hi_fit = theilslopes(stand_hi, x=stand_umin, alpha=0.997)
    umin_co10_corr.append(co10_fit)
    umin_co21_corr.append(co21_fit)
    umin_hi_corr.append(hi_fit)

    # ax.scatter(tab['UMIN'][mask], tab['HI'][mask], label="HI", alpha=0.5)
    # ax.scatter(tab['UMIN'][mask], tab['CO10'][mask], label="CO10", alpha=0.5)
    ax.scatter(stand_umin, stand_hi, label="HI", alpha=0.5,)
    ax.plot(stand_umin, hi_fit[0] * stand_umin + hi_fit[1],)
    ax.scatter(stand_umin, stand_co10, label="CO10", alpha=0.5)
    ax.plot(stand_umin, co10_fit[0] * stand_umin + co10_fit[1])
    ax.scatter(stand_umin, stand_co21, label="CO21", alpha=0.5)
    ax.plot(stand_umin, co21_fit[0] * stand_umin + co21_fit[1])
    ax.grid()

    if i == 0:
        ax.legend(frameon=True)

fig.text(0.5, 0.04, 'Stand. Umin', ha='center', va='center')
fig.text(0.06, 0.5, 'Stand. Intensity',
         ha='center', va='center', rotation='vertical')
fig.savefig(osjoin(fig_path, "umin_intensity_per_Rgal.png"))
fig.savefig(osjoin(fig_path, "umin_intensity_per_Rgal.pdf"))
plt.close()

umin_co10_corr = np.array(umin_co10_corr)
umin_co21_corr = np.array(umin_co21_corr)
umin_hi_corr = np.array(umin_hi_corr)

plt.figure()
(_, caps, _) = plt.errorbar(np.arange(7) + 0.5, umin_hi_corr[:, 0],
                            yerr=[umin_hi_corr[:, 0] - umin_hi_corr[:, 2],
                                  umin_hi_corr[:, 3] - umin_hi_corr[:, 0]],
                            label='HI', drawstyle='steps-mid', capsize=10)
for cap in caps:
    cap.set_markeredgewidth(3)
(_, caps, _) = plt.errorbar(np.arange(7) + 0.5, umin_co10_corr[:, 0],
                            yerr=[umin_co10_corr[:, 0] - umin_co10_corr[:, 2],
                                  umin_co10_corr[:, 3] - umin_co10_corr[:, 0]],
                            label='CO10', drawstyle='steps-mid', capsize=10)
for cap in caps:
    cap.set_markeredgewidth(3)
(_, caps, _) = plt.errorbar(np.arange(7) + 0.5, umin_co21_corr[:, 0],
                            yerr=[umin_co21_corr[:, 0] - umin_co21_corr[:, 2],
                                  umin_co21_corr[:, 3] - umin_co21_corr[:, 0]],
                            label='CO21', drawstyle='steps-mid', capsize=10)
for cap in caps:
    cap.set_markeredgewidth(3)
plt.xlabel("Rgal (kpc)")
plt.ylabel("Correlation")
plt.legend(frameon=True)
plt.grid()
plt.savefig(osjoin(fig_path, "umin_intensity_corr_per_Rgal.png"))
plt.savefig(osjoin(fig_path, "umin_intensity_corr_per_Rgal.pdf"))
plt.close()

# Try masking to keep only the bright CO(2-1)
# mask = tab["CO21"] > 0.05
mask = np.isfinite(tab["CO21"]) & (tab['CO10'] > 0)
# mask = np.isfinite(tab["CO21"]) & (tab["CO21"] > 0.8) & (tab['CO10'] > 0)

gamma = tab['GAMMA'][mask]
gamma[gamma == 0.0] = 1e-5

# Make a corner plot
corner(np.array([tab["CO21"][mask], tab["CO10"][mask], tab["HI"][mask],
                 tab["UMIN"][mask], np.log10(gamma)]).T,
       bins=12, data_kwargs={"alpha": 0.7},
       labels=['CO21', 'CO10', 'HI', 'UMIN', "log10 GAMMA"])
plt.subplots_adjust(hspace=.03, wspace=.03)
plt.savefig(osjoin(fig_path, "umin_intensity_corner.png"))
plt.savefig(osjoin(fig_path, "umin_intensity_corner.pdf"))
plt.close()


oldmask = np.isfinite(oldtab["CO21"]) & (oldtab['CO10'] > 0) & \
    np.isfinite(oldtab['HI'])

corner(np.array([tab["CO21"][oldmask], tab["CO10"][oldmask],
                 tab["HI"][oldmask] / 1000.,
                 tab["UMIN"][oldmask]]).T,
       bins=12, data_kwargs={"alpha": 0.7},
       labels=['CO21', 'CO10', 'HI', 'UMIN'])
plt.subplots_adjust(hspace=.03, wspace=.03)

# Standardize the data and calculate the distance from the one-to-one line
stand_umin = standardize(np.array(tab['UMIN'][mask]))
stand_co21 = standardize(np.array(tab['CO21'][mask]))
stand_co10 = standardize(np.array(tab['CO10'][mask]))
stand_hi = standardize(np.array(tab['HI'][mask]))

# Remove the NGC 604 outliers
# outlier_mask = (stand_umin < 2) & (stand_co21 < 2)
outlier_mask = ~np.logical_and(stand_umin > 4, stand_co21 > 4.5)
outlier_mask = outlier_mask & (stand_co21 > 0)


stand_umin = stand_umin[outlier_mask]
stand_co21 = stand_co21[outlier_mask]
stand_co10 = stand_co10[outlier_mask]
stand_hi = stand_hi[outlier_mask]

# Calculate robust correlations
umin_co10_fit = theilslopes(stand_co10, x=stand_umin, alpha=0.997)
umin_co21_fit = theilslopes(stand_co21, x=stand_umin, alpha=0.997)
umin_hi_fit = theilslopes(stand_hi, x=stand_umin, alpha=0.997)

dists21 = line_distance(stand_umin, stand_co21, *umin_co21_fit[:2])
dists10 = line_distance(stand_umin, stand_co10, *umin_co10_fit[:2])
distshi = line_distance(stand_umin, stand_hi, *umin_hi_fit[:2])

# norm = mpl.colors.Normalize(vmin=-3., vmax=3.)

# beam = Beam(60. * u.arcsec)

# co21_mom0 = Projection.from_hdu(fits.open(iram_co21_data_path("m33.co21_iram.mom0.fits"))[0])
# co21_rms = Projection.from_hdu(fits.open(iram_co21_data_path("m33.rms.masked.fits"))[0])

# smooth_co21 = co21_mom0.convolve_to(beam)
# # Remove regions outside of the original map extent
# smooth_co21[np.isnan(co21_rms)] = np.NaN

# smooth_co21.quicklook()
# smooth_co21.FITSFigure.show_markers(tab['RA'][mask],
#                                     tab['DEC'][mask],
#                                     c=dists21, cmap='seismic',
#                                     norm=norm)
# smooth_co21.FITSFigure.savefig(osjoin(fig_path, "m33_co21_corrdist_umin_vs_co21.png"))
# smooth_co21.FITSFigure.savefig(osjoin(fig_path, "m33_co21_corrdist_umin_vs_co21.pdf"))
# smooth_co21.FITSFigure.close()

# del smooth_co21.FITSFigure
# smooth_co21.quicklook()
# smooth_co21.FITSFigure.show_markers(tab['RA'][mask],
#                                     tab['DEC'][mask],
#                                     c=dists10, cmap='seismic',
#                                     norm=norm)
# smooth_co21.FITSFigure.savefig(osjoin(fig_path, "m33_co21_corrdist_umin_vs_co10.png"))
# smooth_co21.FITSFigure.savefig(osjoin(fig_path, "m33_co21_corrdist_umin_vs_co10.pdf"))
# smooth_co21.FITSFigure.close()

# del smooth_co21.FITSFigure
# smooth_co21.quicklook()
# smooth_co21.FITSFigure.show_markers(tab['RA'],
#                                     tab['DEC'],
#                                     c=tab['UMIN'], cmap='viridis',
#                                     norm=mpl.colors.Normalize(vmin=0, vmax=5))
# smooth_co21.FITSFigure.savefig(osjoin(fig_path, "m33_co21_umin.png"))
# smooth_co21.FITSFigure.savefig(osjoin(fig_path, "m33_co21_umin.pdf"))
# smooth_co21.FITSFigure.close()
