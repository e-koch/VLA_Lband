
'''
Compare the HI properties w/ and w/o CO detected
'''

import numpy as np
import astropy.units as u
import matplotlib.pyplot as p
import astropy.io.fits as fits
from astropy.table import Table
from radio_beam import Beam
from scipy.stats import binned_statistic
from scipy import ndimage as nd
import os
from os.path import exists
from os.path import join as osjoin
from corner import hist2d
from skimage.segmentation import find_boundaries
from astropy.visualization import hist as astro_hist
import seaborn as sb

from paths import (iram_co21_14B088_data_path, fourteenB_HI_data_wGBT_path,
                   fourteenB_wGBT_HI_file_dict,
                   allfigs_path)
from constants import (co21_mass_conversion, hi_mass_conversion,
                       hi_freq, ang_to_phys)
from galaxy_params import gal_feath as gal
from plotting_styles import (onecolumn_figure, default_figure,
                             twocolumn_figure, twocolumn_twopanel_figure)

fig_path = osjoin(allfigs_path(""), "co_vs_hi")
if not exists(fig_path):
    os.mkdir(fig_path)

# Load in radial averages
dr = 100.0 * u.pc
tab_name = "tables/co21_hi_radialprofiles_{}pc.fits".format(int(dr.value))
try:
    tab = Table.read(fourteenB_HI_data_wGBT_path(tab_name))
except OSError:
    raise OSError("Table does not exist in the 14B-088 data path. "
                  "Run co_radial_profile.py first:"
                  " {}".format(tab_name))

rs = u.Quantity(tab["Radius"])
sd = u.Quantity(tab["CO_Sigma"])
sd_sigma = u.Quantity(tab["CO_Sigma_std"])
sd_hi = u.Quantity(tab["HI_Sigma"])
sd_sigma_hi = u.Quantity(tab["HI_Sigma_std"])

# Now plot their ratio against the total gas surface density
gas_ratio = sd.value / sd_hi.value
gas_ratio_sigma = \
    (gas_ratio *
     np.sqrt((sd_sigma / sd)**2 + (sd_sigma_hi / sd_hi)**2)).value
log_gas_ratio_sigma = gas_ratio_sigma / (gas_ratio * np.log(10))
total_sd = sd.value + sd_hi.value
total_sd_sigma = \
    (total_sd *
     np.sqrt((sd_sigma / sd)**2 + (sd_sigma_hi / sd_hi)**2)).value

# Load in the same data created in h2_hi_ratios.py

hi_mom0 = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]
# hi_mom0 = fits.open(fourteenB_HI_file_dict["Moment0"])[0]
beam = Beam.from_fits_header(hi_mom0.header)
# Convert to K km/s
hi_mom0_data = hi_mom0.data * beam.jtok(hi_freq).value / 1000.
hi_mom0_data = hi_mom0_data * u.K * u.km / u.s

hi_lwidth = fits.open(fourteenB_wGBT_HI_file_dict['LWidth'])[0].data
hi_skew = fits.open(fourteenB_wGBT_HI_file_dict['Skewness'])[0].data
hi_kurt = fits.open(fourteenB_wGBT_HI_file_dict['Kurtosis'])[0].data
hi_peaktemp = fits.open(fourteenB_wGBT_HI_file_dict['PeakTemp'])[0].data

co_sig_mask = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_source_mask.fits"))[0]

# good_pts = np.where(np.isfinite(mom0_reproj))
co_mask = co_sig_mask.data.sum(0) >= 2
good_pts = np.where(co_mask)

# Load in the per-pixel column densities
tab = Table.read(fourteenB_HI_data_wGBT_path("tables/column_densities_perpix.fits"))

hi_coldens = tab['Sigma_HI'] * u.solMass / u.pc**2
co_coldens = tab['Sigma_H2'] * u.solMass / u.pc**2
radii_pts = tab['Radius'] * u.kpc
pang_pts = tab['PA'] * u.deg
gas_ratio_pix = tab['Ratio'] * u.dimensionless_unscaled
total_sd_pix = tab['Sigma_Total'] * u.solMass / u.pc**2


# Compare different properties vs. radius and PA
twocolumn_figure()

fig, ax = p.subplots(2, 2, sharex=True)
hist2d(radii_pts.value, np.log10(gas_ratio_pix.value),
       data_kwargs={"alpha": 0.6}, ax=ax[0, 0], bins=10)
# Make radial bins and find the median in each
# rad_bins = np.arange(0, 7.5, 0.5)
# med_ratio, bin_edges, cts = binned_statistic(radii_pts.value,
#                                              np.log10(gas_ratio_pix.value),
#                                              bins=rad_bins,
#                                              statistic=np.median)
# lower_ratio = binned_statistic(radii_pts.value,
#                                np.log10(gas_ratio_pix.value),
#                                bins=rad_bins,
#                                statistic=lambda x: np.percentile(x, 15))[0]
# upper_ratio = binned_statistic(radii_pts.value,
#                                np.log10(gas_ratio_pix.value),
#                                bins=rad_bins,
#                                statistic=lambda x: np.percentile(x, 85))[0]
# bin_cents = (bin_edges[1:] + bin_edges[:-1]) / 2.
# ax[0, 0].errorbar(bin_cents, med_ratio, fmt='o--',
#                   yerr=[med_ratio - lower_ratio, upper_ratio - med_ratio])

# Overplot the 100 pc radial averages
ax[0, 0].errorbar(rs.value, np.log10(gas_ratio),
                  yerr=log_gas_ratio_sigma,
                  alpha=0.6, fmt='D-')

ax[0, 0].set_ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
                    " \Sigma_{\mathrm{HI}}$")
ax[0, 0].grid()

# ax[1].plot(radii_pts.value, total_sd_pix, "ko",
#            ms=3.0, mec=None, alpha=0.8)
hist2d(radii_pts.value, total_sd_pix.value,
       data_kwargs={"alpha": 0.6}, ax=ax[0, 1], bins=12)
ax[0, 1].errorbar(rs.value, total_sd,
                  yerr=total_sd_sigma,
                  alpha=0.6, fmt='D-')

ax[0, 1].set_ylabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
ax[0, 1].grid()

hist2d(radii_pts.value, co_coldens.value,
       data_kwargs={"alpha": 0.6}, ax=ax[1, 0], bins=12)
ax[1, 0].set_ylabel("$\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
ax[1, 0].set_xlabel("Radius (kpc)")
ax[1, 0].grid()
ax[1, 0].errorbar(rs.value, sd.value,
                  yerr=sd_sigma.value,
                  alpha=0.6, fmt='D-')
hist2d(radii_pts.value, hi_coldens.value,
       data_kwargs={"alpha": 0.6}, ax=ax[1, 1], bins=12)
ax[1, 1].set_ylabel("$\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
ax[1, 1].set_xlabel("Radius (kpc)")
ax[1, 1].grid()
ax[1, 1].errorbar(rs.value, sd_hi.value,
                  yerr=sd_sigma_hi.value,
                  alpha=0.6, fmt='D-')

ax[0, 0].set_xlim(0, 7.2)

save_name = "ratio_surfdens_vs_radius_perpix"
p.savefig(osjoin(fig_path, "{0}_{1}pc.pdf".format(save_name,
                                                  int(dr.value))))
p.savefig(osjoin(fig_path, "{0}_{1}pc.png".format(save_name,
                                                  int(dr.value))))
p.close()

# Comparison with radius

twocolumn_figure()

# Only go out to 7 kpc where the map is azimuthally complete
bin_width = 1.2
inners = np.arange(0, 7, bin_width)
outers = np.arange(bin_width, 7.5, bin_width)
fig, axes = p.subplots(3, 2, sharex=True, sharey=True)
p.subplots_adjust(hspace=0.04, wspace=0.04)

for i, (ax, lower, upper) in enumerate(zip(axes.flatten(), inners,
                                           outers)):
    bin_pts = np.logical_and(radii_pts.value >= lower,
                             radii_pts.value < upper)
    hist2d(total_sd_pix.value,
           np.log10(gas_ratio_pix.value),
           data_kwargs={"alpha": 0.3},
           ax=ax)
    ax.plot(total_sd_pix[bin_pts], np.log10(gas_ratio_pix[bin_pts]), "D",
            label="{0}-{1}".format(lower, upper),
            ms=3.0, mec=None, alpha=0.3)
    ax.grid()
    ax.text(40, -0.8, "{0} - {1} kpc".format(lower, upper),
            bbox={"boxstyle": "square", "facecolor": "w"},
            fontsize=11)

ax.set_xlim([0, 75])
ax.set_ylim([-1.5, 0.8])

fig.text(0.5, 0.04, "$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)",
         ha='center', va='center',)
fig.text(0.06, 0.5, "log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} / \Sigma_{\mathrm{HI}}$",
         ha='center', va='center', rotation='vertical',)

save_name = "ratio_totalsigma_radialbins_perpix"
fig.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
fig.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()


# Nothing interesting vs. position angle

# fig, ax = p.subplots(1, 2)
# # hist2d(pang_pts.value, np.log10(gas_ratio_pix.value),
# #        data_kwargs={"alpha": 0.6}, ax=ax[0], bins=10)
# ax[0].plot(pang_pts.value, np.log10(gas_ratio_pix.value), "ko",
#            ms=3.0, mec=None, alpha=0.8)
# ax[0].set_ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
#                  " \Sigma_{\mathrm{HI}}$")
# ax[0].set_xlabel("Position Angle (deg)")
# ax[1].plot(pang_pts.value, total_sd_pix, "ko",
#            ms=3.0, mec=None, alpha=0.8)
# ax[1].set_ylabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
# ax[1].set_xlabel("Position Angle (deg)")

# Something interesting of the above plot is the difference in HI
# properties between the HI azimuthal average and the HI where CO
# is detected. Compare the populations of HI w/ CO and all HI

inc = np.cos(gal.inclination)
all_hi_pts = np.isfinite(hi_mom0_data)
radii = gal.radius(header=hi_mom0.header).to(u.kpc)
hi_coldens_all = hi_mom0_data[all_hi_pts] * hi_mass_conversion * inc
radii_pts_all = radii[all_hi_pts]

onecolumn_figure()

hist2d(radii_pts_all.value, hi_coldens_all.value,
       data_kwargs={"alpha": 0.6}, bins=30)
p.plot(radii_pts.value, hi_coldens.value, 'D', alpha=0.1)

p.ylabel("$\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel("Radius (kpc)")
p.tight_layout()

save_name = "sigma_HI_perpix_wCO_and_all"
p.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
p.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()

# Calculate the fraction of LOS w/ CO as function of HI column density
bins = np.linspace(0, hi_coldens_all.max())
bin_cents = (bins[:-1] + bins[1:]) / 2.
co_fraction = []
for low, high in zip(bins[:-1], bins[1:]):
    num_co = np.logical_and(hi_coldens.value > low.value,
                            hi_coldens.value <= high.value)
    num_total = np.logical_and(hi_coldens_all > low, hi_coldens_all <= high)

    co_fraction.append(num_co.sum() / float(num_total.sum()))

co_fraction = np.array(co_fraction)

# Now do the same for the peak HI temperature
bins_peak = np.linspace(hi_peaktemp[all_hi_pts].min(),
                        hi_peaktemp[all_hi_pts].max())
bin_cents_peak = (bins_peak[:-1] + bins_peak[1:]) / 2.

co_fraction_peak = []
for low, high in zip(bins_peak[:-1], bins_peak[1:]):
    num_co = np.logical_and(hi_peaktemp[good_pts] > low,
                            hi_peaktemp[good_pts] <= high)
    num_total = np.logical_and(hi_peaktemp[all_hi_pts] > low,
                               hi_peaktemp[all_hi_pts] <= high)

    co_fraction_peak.append(num_co.sum() / float(num_total.sum()))

co_fraction_peak = np.array(co_fraction_peak)

onecolumn_figure()

ax = p.subplot(111)
pl1 = ax.plot(bin_cents_peak, co_fraction_peak, '--',
              color=sb.color_palette()[1], drawstyle='steps-mid',
              label=r"T$_{\rm HI}$")
ax.grid()
ax.set_xlabel(r"T$_{\rm HI}$ (K)")
ax.set_ylabel("CO Detection Fraction")

ax2 = ax.twiny()

pl2 = ax2.plot(bin_cents, co_fraction, drawstyle='steps-mid',
               label=r'$\Sigma_{\mathrm{HI}}$')
ax2.set_xlabel("$\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")

pls = [pl2[0], pl1[0]]
labs = [r'$\Sigma_{\mathrm{HI}}$', r"T$_{\rm HI}$"]
ax2.legend(pls, labs, frameon=True, loc='upper left')

p.tight_layout()

save_name = "sigma_peak_HI_perpix_CO_detection_fraction"
p.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
p.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()

# Compare the detection fraction as fcn of radius
twocolumn_figure()

fig, axes = p.subplots(3, 2, sharex=True, sharey=True)
p.subplots_adjust(hspace=0.04, wspace=0.04)

co_fractions = []
co_peak_fractions = []

for i, (ax, lower, upper) in enumerate(zip(axes.flatten(), inners,
                                           outers)):
    bin_pts = np.logical_and(radii_pts.value >= lower,
                             radii_pts.value < upper)
    all_rad_pts = np.logical_and(radii_pts_all.value >= lower,
                                 radii_pts_all.value < upper)
    all_rad_pts_good = np.logical_and(radii[good_pts].value >= lower,
                                      radii[good_pts].value < upper)

    co_fraction_rad = []
    for low, high in zip(bins[:-1], bins[1:]):
        num_co = np.logical_and(hi_coldens[bin_pts].value > low.value,
                                hi_coldens[bin_pts].value <= high.value)
        num_total = np.logical_and(hi_coldens_all[all_rad_pts] > low,
                                   hi_coldens_all[all_rad_pts] <= high)

        if num_total.sum() == 0:
            co_fraction_rad.append(np.NaN)
        else:
            co_fraction_rad.append(num_co.sum() / float(num_total.sum()))

    co_fraction_rad = np.array(co_fraction_rad)
    co_fractions.append(co_fraction_rad)

    co_fraction_peak_rad = []
    for low, high in zip(bins_peak[:-1], bins_peak[1:]):
        num_co = np.logical_and(hi_peaktemp[good_pts][all_rad_pts_good] > low,
                                hi_peaktemp[good_pts][all_rad_pts_good] <= high)
        num_total = np.logical_and(hi_peaktemp[all_hi_pts][all_rad_pts] > low,
                                   hi_peaktemp[all_hi_pts][all_rad_pts] <= high)

        if num_total.sum() == 0:
            co_fraction_peak_rad.append(np.NaN)
        else:
            co_fraction_peak_rad.append(num_co.sum() / float(num_total.sum()))

    co_fraction_peak_rad = np.array(co_fraction_peak_rad)

    co_peak_fractions.append(co_fraction_peak_rad)

    pl1 = ax.plot(bin_cents_peak, co_fraction_peak_rad, '--',
                  color=sb.color_palette()[1], drawstyle='steps-mid',
                  label=r"T$_{\rm HI}$")
    ax.grid()

    ax2 = ax.twiny()

    pl2 = ax2.plot(bin_cents, co_fraction_rad, drawstyle='steps-mid',
                   label=r'$\Sigma_{\mathrm{HI}}$')

    if i == 2:
        pls = [pl2[0], pl1[0]]
        labs = [r'$\Sigma_{\mathrm{HI}}$', r"T$_{\rm HI}$"]
        ax2.legend(pls, labs, frameon=True, loc='center left')

    ax.text(5, 0.8, "{0} - {1} kpc".format(lower, upper),
            bbox={"boxstyle": "square", "facecolor": "w"},
            fontsize=11)

    if i != 0 and i != 1:
        ax2.tick_params(labeltop='off')

fig.text(0.5, 0.04, r"T$_{\rm HI}$ (K)",
         ha='center', va='center',)
fig.text(0.5, 0.96, r"$\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)",
         ha='center', va='center',)
fig.text(0.06, 0.5, "CO Detection Fraction",
         ha='center', va='center', rotation='vertical',)

save_name = "sigma_peak_HI_perpix_CO_detection_fraction_radialbins"
fig.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
fig.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()

# What does radius vs Sigma_HI look like around the edge of the clouds?

# hist2d(radii_pts.value, hi_coldens.value,
#        bins=10, data_kwargs={"alpha": 0.6})

edge_mask = find_boundaries(co_mask, connectivity=2, mode='outer')
co_dist_mask = nd.distance_transform_edt(~edge_mask)
co_dist_mask[co_mask] = - co_dist_mask[co_mask]

# # Include all pixels within 1 beam width
# boundary_mask = np.logical_and(co_dist_mask > 0, co_dist_mask <= 14)
# boundary_mask = np.logical_and(all_hi_pts, boundary_mask)

# hi_coldens_co_bound = \
#     hi_mom0_data[boundary_mask] * hi_mass_conversion * inc

# hist2d(radii[boundary_mask].value, hi_coldens_co_bound.value,
#        bins=10, data_kwargs={"alpha": 0.6})

# # boundary_mask = np.logical_and(co_dist_mask > 6.5, co_dist_mask < 6.5 * 2)
# boundary_mask = (co_dist_mask > 14) & (co_dist_mask < 28)
# boundary_mask = np.logical_and(all_hi_pts, boundary_mask)

# hi_coldens_co_bound = \
#     hi_mom0_data[boundary_mask] * hi_mass_conversion * inc

# hist2d(radii[boundary_mask].value, hi_coldens_co_bound.value,
#        bins=10, data_kwargs={"alpha": 0.6})


# HI column density vs. distance from CO detection
phys_scale = ang_to_phys(hi_mom0.header['CDELT2'] * u.deg)

co_dists = co_dist_mask[all_hi_pts] * phys_scale.value / 1000.

max_rad = 0.4

hist2d(co_dists[co_dists < max_rad], hi_coldens_all[co_dists < max_rad].value,
       bins=30, data_kwargs={"alpha": 0.4})

# p.scatter((co_dist_mask - co_inverse_dist_mask)[all_hi_pts],
#             hi_coldens_all.value,
#             bins=30, data_kwargs={"alpha": 0.4})

# Make radial bins and find the median in each
rad_bins = np.arange(-0.2, max_rad, 0.1)
med_ratio, bin_edges, cts = binned_statistic(co_dists,
                                             hi_coldens_all.value,
                                             bins=rad_bins,
                                             statistic=np.median)
lower_ratio = binned_statistic(co_dists,
                               hi_coldens_all.value,
                               bins=rad_bins,
                               statistic=lambda x: np.percentile(x, 15))[0]
upper_ratio = binned_statistic(co_dists,
                               hi_coldens_all.value,
                               bins=rad_bins,
                               statistic=lambda x: np.percentile(x, 85))[0]
bin_cents = (bin_edges[1:] + bin_edges[:-1]) / 2.
p.errorbar(bin_cents, med_ratio, fmt='o--',
           yerr=[med_ratio - lower_ratio, upper_ratio - med_ratio])

p.axvline(0, color=sb.color_palette()[1],
          linewidth=3, linestyle='--', alpha=0.8)

p.xlabel("Distance from CO detection (kpc)")
p.ylabel("$\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")

p.grid()
p.tight_layout()

save_name = "sigma_HI_vs_distance_from_CO"
p.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
p.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()

# Binned by radius
# bin_width = 1
# inners = np.arange(0, 8, bin_width)
# outers = np.arange(bin_width, 8.5, bin_width)
# fig, axes = p.subplots(2, 4, sharex=True, sharey=True)
# p.subplots_adjust(hspace=0.04, wspace=0.04)

# twocolumn_figure()

# for i, (ax, lower, upper) in enumerate(zip(axes.flatten(), inners, outers)):
#     bin_pts = np.logical_and(radii[all_hi_pts].value > lower,
#                              radii[all_hi_pts].value <= upper)
#     hist2d(co_dists[bin_pts],
#            hi_coldens_all[bin_pts].value,
#            data_kwargs={"alpha": 0.6},
#            ax=ax)
#     ax.grid()
#     ax.text(40, -0.8, "{0} - {1} kpc".format(int(lower), int(upper)),
#             bbox={"boxstyle": "square", "facecolor": "w"})
#     # ax.set_xlim([0, 60])
#     # ax.set_ylim([-2.1, 0.7])

#     if i == 0 or i == 2:
#         ax.set_ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
#                       " \Sigma_{\mathrm{HI}}$")
#     if i == 2 or i == 3:
#         ax.set_xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")

# Are there other HI line properties that correlate with CO detections?

# L Width (from 2nd moment)
# _ = astro_hist(hi_lwidth[np.isfinite(hi_lwidth)] / 1000., bins='scott',
#                alpha=0.7, normed=True)
# _ = astro_hist(hi_lwidth[good_pts] / 1000., bins='scott', alpha=0.7,
#                normed=True)

# The distributions are basically indistinguishable. Not surprising given how
# inaccurate the HI line width is

# Skewness
# _ = astro_hist(hi_skew[np.isfinite(hi_skew)], bins='scott',
#                alpha=0.7, normed=True)
# _ = astro_hist(hi_skew[good_pts], bins='scott', alpha=0.7,
#                normed=True)

# Same. Slight preference for small negative skew, but that's probably
# unassociated with the CO

# Kurtosis
# _ = astro_hist(hi_kurt[np.isfinite(hi_kurt)], bins='scott',
#                alpha=0.7, normed=True)
# _ = astro_hist(hi_kurt[good_pts], bins='scott', alpha=0.7,
#                normed=True)

# Distinct preference for positive kurtosis. Now check if this is driven
# by kurtosis - HI peak correlation

# Peak temperature
# _ = astro_hist(hi_peaktemp[np.isfinite(hi_peaktemp)], bins='scott',
#                alpha=0.7, normed=True)
# _ = astro_hist(hi_peaktemp[good_pts], bins='scott', alpha=0.7,
#                normed=True)

# Most definitely! So the kurtosis population difference isn't saying anything
# about the line shape around the CO component.

default_figure()
