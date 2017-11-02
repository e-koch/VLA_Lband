
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

from paths import (iram_co21_14B088_data_path, fourteenB_HI_data_wGBT_path,
                   fourteenB_wGBT_HI_file_dict,
                   allfigs_path)
from constants import (co21_mass_conversion, hi_mass_conversion,
                       hi_freq, ang_to_phys)
from galaxy_params import gal
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

mom0_reproj = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits"))[0]
mom0_reproj = (mom0_reproj.data / 1000.) * u.K * u.km / u.s

good_pts = np.where(np.isfinite(mom0_reproj))

# Make a radius array
radii = gal.radius(header=hi_mom0.header).to(u.kpc)
radii_pts = radii[good_pts]

# And the position angles
pang = gal.position_angles(header=hi_mom0.header).to(u.deg)
pang_pts = pang[good_pts]

skycoord_grid = gal.skycoord_grid(header=hi_mom0.header)
skycoord_pts = skycoord_grid[good_pts]

# Correct for the disk inclincation
inc = np.cos(gal.inclination)

# 30 m beam efficiency
beam_eff = 0.75

# Convert the integrated intensities to surface densities.
# hi_coldens = hi_mom0_reproj[good_pts] * hi_mass_conversion * inc
hi_coldens = hi_mom0_data[good_pts] * hi_mass_conversion * inc

# co_coldens = mom0[good_pts] * co21_mass_conversion * inc / beam_eff
co_coldens = mom0_reproj[good_pts] * co21_mass_conversion * inc / beam_eff

# Remove any NaNs in either
nans = np.logical_or(np.isnan(hi_coldens), np.isnan(co_coldens))

hi_coldens = hi_coldens[~nans]
co_coldens = co_coldens[~nans]
radii_pts = radii_pts[~nans]
pang_pts = pang_pts[~nans]
skycoord_pts = skycoord_pts[~nans]
ypts = good_pts[0][~nans]
xpts = good_pts[1][~nans]

gas_ratio_pix = co_coldens / hi_coldens
total_sd_pix = co_coldens + hi_coldens


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
    ax.set_xlim([0, 60])
    ax.set_ylim([-2.1, 0.7])

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

all_hi_pts = np.isfinite(hi_mom0_data)
hi_coldens_all = hi_mom0_data[all_hi_pts] * hi_mass_conversion * inc
radii_pts_all = radii[all_hi_pts]

onecolumn_figure()

hist2d(radii_pts_all.value, hi_coldens_all.value,
       data_kwargs={"alpha": 0.6}, bins=30)
p.scatter(radii_pts.value, hi_coldens.value, alpha=0.1)

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
    num_co = np.logical_and(hi_coldens > low, hi_coldens <= high)
    num_total = np.logical_and(hi_coldens_all > low, hi_coldens_all <= high)

    co_fraction.append(num_co.sum() / float(num_total.sum()))

co_fraction = np.array(co_fraction)
p.plot(bin_cents, co_fraction)
p.grid()
p.xlabel("$\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("CO Detection Fraction")
p.tight_layout()

save_name = "sigma_HI_perpix_CO_detection_fraction"
p.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
p.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()

# What does radius vs Sigma_HI look like around the edge of the clouds?

# hist2d(radii_pts.value, hi_coldens.value,
#        bins=10, data_kwargs={"alpha": 0.6})

co_mask = np.isfinite(mom0_reproj)
co_dist_mask = nd.distance_transform_edt(~co_mask)
# co_inverse_dist_mask = nd.distance_transform_edt(co_mask)

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

hist2d(co_dists, hi_coldens_all.value,
       bins=30, data_kwargs={"alpha": 0.4})

# p.scatter((co_dist_mask - co_inverse_dist_mask)[all_hi_pts],
#             hi_coldens_all.value,
#             bins=30, data_kwargs={"alpha": 0.4})

# Make radial bins and find the median in each
rad_bins = np.arange(0, 4, 0.2)
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

default_figure()
