
import matplotlib.pyplot as p
from astropy.table import Table
import os
from os.path import exists
from os.path import join as osjoin
from corner import hist2d
import seaborn as sb
from scipy.stats import binned_statistic
import astropy.units as u

from cube_analysis.h2_models import (krumholz2013_ratio_model,
                                     krumholz2013_sigmaHI,
                                     optimize_clump_factors)

from paths import (fourteenB_HI_data_wGBT_path,
                   allfigs_path)
from plotting_styles import (onecolumn_figure, default_figure,
                             twocolumn_figure, twocolumn_twopanel_figure)

fig_path = osjoin(allfigs_path(""), "co_vs_hi/h2_formation_models")
if not exists(fig_path):
    os.mkdir(fig_path)

# Bin radius should match whatever was used in co_radial_profile.py
dr = 100.0 * u.pc

# Read in the radial profiles table
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

# Overplot the Krumholz model with a few different clumping factors.
# Theoretically, c -> 1 at a resolution of 100 pc. but I'm finding a better
# match when c=4-6. The model is supposed to take metallicity into account,
# but maybe the gradient is causing some issues? Schruba+11 finds c~2 for
# their entire sample, with a lot of scatter
sds = np.arange(1, 40, 0.2)

onecolumn_figure(font_scale=1.0)
# p.semilogy(total_sd, gas_ratio, 'bD')
p.errorbar(total_sd, np.log10(gas_ratio), yerr=log_gas_ratio_sigma,
           xerr=total_sd_sigma, alpha=0.6, fmt='D')
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=1, Z=1.0)), "--",
       label="c=1, Z=1.0")
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=2, Z=0.5)), "-.",
       label="c=2, Z=0.5")
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=0.5)), ":",
       label="c=3, Z=0.5")
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=1.0)), "-",
       label="c=3, Z=1.0")
p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
         " \Sigma_{\mathrm{HI}}$")
p.xlim([2, 22])
p.ylim([-4, 1])
p.legend(loc='lower right', frameon=True)
p.grid()
p.tight_layout()

save_name = "ratio_totalsigma_w_krumholzmodel_dr"
p.savefig(osjoin(fig_path, "{0}_{1}pc.pdf".format(save_name,
                                                  int(dr.value))))
p.savefig(osjoin(fig_path, "{0}_{1}pc.png".format(save_name,
                                                  int(dr.value))))
p.close()

# Gratier+16 find evidence for a dark CO component, at about ~5 Msol/pc^2.
# Let's add this in, assuming the dark component is *only* in the CO and
# not due to optically thick HI (some portion probably is).
sd_dark = sd + 5 * u.solMass / u.pc**2
sd_dark_sigma = (sd_dark * sd_sigma) / sd
gas_ratio_dark = sd_dark.value / sd_hi.value
gas_ratio_dark_sigma = \
    (gas_ratio_dark *
     np.sqrt((sd_dark_sigma / sd_dark)**2 +
             (sd_sigma_hi / sd_hi)**2)).value
log_gas_ratio_dark_sigma = gas_ratio_dark_sigma / \
    (gas_ratio_dark * np.log(10))

total_sd_plus_dark = sd_dark.value + sd_hi.value
total_sd_plus_dark_sigma = \
    (total_sd_plus_dark *
     np.sqrt((sd_dark_sigma / sd_dark)**2 +
             (sd_sigma_hi / sd_hi)**2)).value

p.plot(total_sd, np.log10(gas_ratio), 'D',
       alpha=0.6, label=r"H$_2$ + HI")
# p.errorbar(total_sd, np.log10(gas_ratio), yerr=log_gas_ratio_sigma,
#            xerr=total_sd_sigma, alpha=0.6, fmt='D',
#            label=r"H$_2$ + HI")
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=1, Z=1.0)), "--",
       label="c=1, Z=1.0")
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=2, Z=0.5)), "-.",
       label="c=2, Z=0.5")
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=0.5)), ":",
       label="c=3, Z=0.5")
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=1.0)), "-",
       label="c=3, Z=1.0")
p.plot(total_sd_plus_dark, np.log10(gas_ratio_dark), 'o',
       alpha=0.6, label=r"H$_2$ + HI + CO-dark H$_2$")
# p.errorbar(total_sd_plus_dark, np.log10(gas_ratio_dark),
#            yerr=log_gas_ratio_dark_sigma,
#            xerr=total_sd_plus_dark_sigma, alpha=0.6, fmt='o',
#            label=r"H$_2$ + HI + CO-dark H$_2$")
p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
         " \Sigma_{\mathrm{HI}}$")
p.xlim([2, 25])
p.ylim([-4, 1])
p.legend(loc='lower right', frameon=True)
p.grid()
p.tight_layout()

save_name = "ratio_totalsigma_dark_w_krumholzmodel_dr"
p.savefig(osjoin(fig_path, "{0}_{1}pc.pdf".format(save_name,
                                                  int(dr.value))))
p.savefig(osjoin(fig_path, "{0}_{1}pc.png".format(save_name,
                                                  int(dr.value))))
p.close()

# But M33 has a known metallicity gradient, so we can do a bit better
# Clumping factors should converge to 1 on 100 pc, based on the Krumholz
# model. This isn't happening here, so let's what c needs to be for the
# curve to intersect with each point we have.

# Metallicity of 0.5
clump_constz = optimize_clump_factors(total_sd, gas_ratio, Z=0.5,
                                      c_init=3.5)

# Metallicity Gradient from Roso & Simon (2005)
def ros_sim_metallicity(radius):
    return 10 ** (8.36 - 0.027 * radius - 8.8)

clump_rossim = optimize_clump_factors(total_sd, gas_ratio,
                                      Z=ros_sim_metallicity(rs.value),
                                      c_init=7.0)

# And from Bresolin 2011, Equation 1
def bresolin_metallicity(radius):
    return 10 ** (8.48 - 0.042 * radius - 8.8)
    # return 10 ** (8.82 - 0.03 * radius - 8.8)

clump_bresolin = optimize_clump_factors(total_sd, gas_ratio,
                                        Z=bresolin_metallicity(rs.value),
                                        c_init=7.)

onecolumn_figure(font_scale=1.0)

p.plot(rs.value[:-1], clump_constz[:-1], 'D-', label="Z=0.5")
p.plot(rs.value[:-1], clump_rossim[:-1], 'o--',
       label="Rosolowsky \& Simon (2008)")
p.plot(rs.value[:-1], clump_bresolin[:-1], 's-.', label="Bresolin (2011)")
p.legend(loc='best', frameon=True)
p.ylim([2, 8])
p.grid()
p.ylabel("Clumping Factor")
p.xlabel("Radius (kpc)")
p.tight_layout()

save_name = "clumpfactor_krumholzmodel_dr"
p.savefig(osjoin(fig_path, "{0}_{1}pc.pdf".format(save_name,
                                                  int(dr.value))))
p.savefig(osjoin(fig_path, "{0}_{1}pc.png".format(save_name,
                                                  int(dr.value))))
p.close()

# What are the properties like on a per pixel basis?

# Load in the per-pixel column densities
tab = Table.read(fourteenB_HI_data_wGBT_path("tables/column_densities_perpix.fits"))

hi_coldens = tab['Sigma_HI'] * u.solMass / u.pc**2
co_coldens = tab['Sigma_H2'] * u.solMass / u.pc**2
radii_pts = tab['Radius'] * u.kpc
pang_pts = tab['PA'] * u.deg
gas_ratio_pix = tab['Ratio'] * u.dimensionless_unscaled
total_sd_pix = tab['Sigma_Total'] * u.solMass / u.pc**2

sds = np.arange(0.1, 75, 0.2)

onecolumn_figure()

hist2d(total_sd_pix.value, np.log10(gas_ratio_pix.value),
       data_kwargs={"alpha": 0.3})
p.xlim([0, 75])
p.ylim([-1.4, 0.8])
p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
         " \Sigma_{\mathrm{HI}}$")

p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=1, Z=0.5)), "-",
       label="c=1, Z=0.5", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=1, Z=1.0)), "--",
       label="c=1, Z=1.0", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=2, Z=0.5)), "-.",
       label="c=2, Z=0.5", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=0.5)), ":",
       label="c=3, Z=0.5", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=1.0)), "-",
       label="c=3, Z=1.0", linewidth=2, alpha=0.95)

p.legend(loc='lower right', frameon=True)
p.grid()

p.tight_layout()

save_name = "ratio_totalsigma_w_krumholzmodel_perpix"
p.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
p.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()

# Overplot the radial averages

hist2d(total_sd_pix.value, np.log10(gas_ratio_pix.value),
       data_kwargs={"alpha": 0.3})
p.xlim([0, 75])
p.ylim([-1.8, 0.8])
p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
         " \Sigma_{\mathrm{HI}}$")

p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=1, Z=0.5)), "-",
       label="c=1, Z=0.5", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=1, Z=1.0)), "--",
       label="c=1, Z=1.0", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=2, Z=0.5)), "-.",
       label="c=2, Z=0.5", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=0.5)), ":",
       label="c=3, Z=0.5", linewidth=2, alpha=0.95)
p.plot(sds, np.log10(krumholz2013_ratio_model(sds, c=3, Z=1.0)), "-",
       label="c=3, Z=1.0", linewidth=2, alpha=0.95)

p.plot(total_sd, np.log10(gas_ratio), 'D', markeredgecolor='k',
       markeredgewidth=0.25,)

p.legend(loc='lower right', frameon=True)
p.grid()

p.tight_layout()

save_name = "ratio_totalsigma_w_krumholzmodel_perpix_w_radavg"
p.savefig(osjoin(fig_path, "{0}_{1}pc.pdf".format(save_name,
                                                  int(dr.value))))
p.savefig(osjoin(fig_path, "{0}_{1}pc.png".format(save_name,
                                                  int(dr.value))))
p.close()

# Sigma HI vs total
cpal = sb.color_palette()

hist2d(total_sd_pix.value, hi_coldens.value,
       data_kwargs={"alpha": 0.3})
p.xlim([0, 75])
p.ylim([-3, 27])
p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("$\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")

# p.axhline(krumholz_maxhi_sigma(c=1, Z=1.0).value, linestyle="--",
#           label="c=1, Z=1.0", linewidth=2, alpha=0.95,
#           color=cpal[1])
# p.axhline(krumholz_maxhi_sigma(c=2, Z=0.5).value, linestyle="-.",
#           label="c=2, Z=0.5", linewidth=2, alpha=0.95,
#           color=cpal[2])
# p.axhline(krumholz_maxhi_sigma(c=3, Z=0.5).value, linestyle=":",
#           label="c=3, Z=0.5", linewidth=2, alpha=0.95,
#           color=cpal[3])
# p.axhline(krumholz_maxhi_sigma(c=3, Z=1.0).value, linestyle="-",
#           label="c=3, Z=1.0", linewidth=2, alpha=0.95,
#           color=cpal[4])

p.plot(sds, krumholz2013_sigmaHI(sds, c=1, Z=1.0), "--",
       label="c=1, Z=1.0", linewidth=2, alpha=0.95,
       color=cpal[1])
p.plot(sds, krumholz2013_sigmaHI(sds, c=2, Z=0.5), "-.",
       label="c=2, Z=0.5", linewidth=2, alpha=0.95,
       color=cpal[2])
p.plot(sds, krumholz2013_sigmaHI(sds, c=3, Z=0.5), ":",
       label="c=3, Z=0.5", linewidth=2, alpha=0.95,
       color=cpal[3])
p.plot(sds, krumholz2013_sigmaHI(sds, c=3, Z=1.0), "-",
       label="c=3, Z=1.0", linewidth=2, alpha=0.95,
       color=cpal[4])

p.plot(total_sd, sd_hi, 'D', markeredgecolor='k',
       markeredgewidth=0.25, color=cpal[0])

p.plot([0, 27], [0, 27], '-', linewidth=4, alpha=0.6, color=cpal[5])

p.legend(loc='lower right', frameon=True, ncol=2)
p.grid()

p.tight_layout()

save_name = "sigma_hi_vs_total_w_krumholzmodel_perpix_w_radavg"
p.savefig(osjoin(fig_path, "{0}_{1}pc.pdf".format(save_name,
                                                  int(dr.value))))
p.savefig(osjoin(fig_path, "{0}_{1}pc.png".format(save_name,
                                                  int(dr.value))))
p.close()

# Can the metallicity gradient account for the variation with radius?
# Solve for clumping factor for every pixel.
# If the metallicity accounts for the variation, the clumping factors
# should be constant with radius, on average.

clump_per_pix_rossim = \
    optimize_clump_factors(total_sd_pix.value,
                           gas_ratio_pix.value,
                           Z=ros_sim_metallicity(radii_pts.value),
                           c_init=10.)

clump_per_pix_bresolin = \
    optimize_clump_factors(total_sd_pix.value,
                           gas_ratio_pix.value,
                           Z=bresolin_metallicity(radii_pts.value),
                           c_init=10.)

clump_per_pix_constz = \
    optimize_clump_factors(total_sd_pix.value,
                           gas_ratio_pix.value,
                           Z=0.5,
                           c_init=10.)

# Note that plotting the scatter of all clump factors does not show the
# population well, but it REALLY highlights those low SB outlier points!

# Make radial bins and find the median in each
rad_bins = np.arange(0, 7.0, 0.5)
med_ratio, bin_edges, cts = binned_statistic(radii_pts.value,
                                             clump_per_pix_rossim,
                                             bins=rad_bins,
                                             statistic=np.median)
lower_ratio = binned_statistic(radii_pts.value,
                               clump_per_pix_rossim,
                               bins=rad_bins,
                               statistic=lambda x: np.percentile(x, 15))[0]
upper_ratio = binned_statistic(radii_pts.value,
                               clump_per_pix_rossim,
                               bins=rad_bins,
                               statistic=lambda x: np.percentile(x, 85))[0]
bin_cents = (bin_edges[1:] + bin_edges[:-1]) / 2.

med_ratio_bres, bin_edges_bres, cts_bres = \
    binned_statistic(radii_pts.value, clump_per_pix_bresolin,
                     bins=rad_bins,
                     statistic=np.median)
lower_ratio_bres = binned_statistic(radii_pts.value,
                                    clump_per_pix_bresolin,
                                    bins=rad_bins,
                                    statistic=lambda x: np.percentile(x, 15))[0]
upper_ratio_bres = binned_statistic(radii_pts.value,
                                    clump_per_pix_bresolin,
                                    bins=rad_bins,
                                    statistic=lambda x: np.percentile(x, 85))[0]

med_ratio_const, bin_edges_const, cts_const = \
    binned_statistic(radii_pts.value, clump_per_pix_constz,
                     bins=rad_bins,
                     statistic=np.median)
lower_ratio_const = binned_statistic(radii_pts.value,
                                     clump_per_pix_constz,
                                     bins=rad_bins,
                                     statistic=lambda x: np.percentile(x, 15))[0]
upper_ratio_const = binned_statistic(radii_pts.value,
                                     clump_per_pix_constz,
                                     bins=rad_bins,
                                     statistic=lambda x: np.percentile(x, 85))[0]

onecolumn_figure()

p.errorbar(bin_cents, med_ratio_const, fmt='D-',
           yerr=[med_ratio_const - lower_ratio_const,
                 upper_ratio_const - med_ratio_const],
           label='Z=0.5')
p.errorbar(bin_cents, med_ratio, fmt='o--',
           yerr=[med_ratio - lower_ratio, upper_ratio - med_ratio],
           label='Rosolowsky \& Simon (2008)')
p.errorbar(bin_cents, med_ratio_bres, fmt='s-.',
           yerr=[med_ratio_bres - lower_ratio_bres,
                 upper_ratio_bres - med_ratio_bres],
           label='Bresolin (2011)')

p.grid()
p.legend(frameon=True, loc='upper right')

p.ylabel("Clumping Factor")
p.xlabel("Radius (kpc)")

p.ylim([0.8, 6.7])

p.tight_layout()

save_name = "clumpfactor_krumholzmodel_perpix_500pc"
p.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
p.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()


# Compare the ratio of the averages. If the gradient takes into account
# the clump factor changes, it is accounting for the differences.
p.plot(bin_cents, med_ratio / med_ratio_const, 'D-',
       label='Ros. + Simon / Const Z.')
p.plot(bin_cents, med_ratio_bres / med_ratio_const, 'o--',
       label='Bresolin / Const Z.')

p.legend(frameon=True, loc='upper left')

p.ylabel("Clumping Factor Ratio")
p.xlabel("Radius (kpc)")

p.grid()

p.tight_layout()

save_name = "clumpfactor_ratio_krumholzmodel_perpix_500pc"
p.savefig(osjoin(fig_path, "{0}.pdf".format(save_name)))
p.savefig(osjoin(fig_path, "{0}.png".format(save_name)))
p.close()

default_figure()
