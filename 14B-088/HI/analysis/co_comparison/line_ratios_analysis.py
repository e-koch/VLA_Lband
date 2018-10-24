
'''
Looking for HI-CO relationships from the fitting in co_hi_list_fitting.py
'''

import os
from os.path import join as osjoin
from spectral_cube import Projection
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from corner import hist2d
import seaborn as sb
import emcee
from scipy.stats import binned_statistic
from scipy import stats

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, iram_co21_14B088_data_path)
from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_twopanel_figure,
                             twocolumn_figure)
from constants import co21_mass_conversion, hi_mass_conversion
from galaxy_params import gal_feath as gal

fig_path = allfigs_path("co_vs_hi")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)


inc = np.cos(gal.inclination)

# Start with checking the column densities from the fitting.

tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_co_gaussfit_column_densities_perpix.fits"))

# Don't consider the "bad fits" that are probably due to multiple components
good_pts = np.logical_and(~tab['multicomp_flag_HI'],
                          ~tab['multicomp_flag_CO'])

# Our flagging didn't pick up some bad HI fits that are supremely narrow
# Inspecting these profiles by-eye doesn't clearly show why some failed, but
# at least half of the ~10 points have multi-components in CO
# Here's the code for inspecting the outliers
# for y, x in zip(tab['ypts'][np.where(tab['sigma_HI'][good_pts] < 3800)],
#                 tab['xpts'][np.where(tab['sigma_HI'][good_pts] < 3800)]):
#     ax = p.subplot(111)
#     ax.plot(co_cube.spectral_axis.value, co_cube[:, y, x].value, drawstyle='steps-mid')
#     ax2 = ax.twinx()
#     ax2.plot(cube.spectral_axis.value, cube[:, y, x].value, drawstyle='steps-mid')
#     p.draw()
#     raw_input("{0}, {1}?".format(y,x))
#     p.clf()

good_pts = np.logical_and(good_pts,
                          tab["sigma_HI"] > 3800)

# Minimum CO line width of one channel.
good_pts = np.logical_and(good_pts,
                          tab["sigma_CO"] >= 2600)


# Show that the Gaussian model is a decent representation
# twocolumn_twopanel_figure()
cpal = sb.color_palette()

ax1 = plt.subplot(121)

hist2d(tab["coldens_CO_gauss"][good_pts] / (inc * co21_mass_conversion).value,
       tab["coldens_CO_FWHM"][good_pts] / (inc * co21_mass_conversion).value,
       data_kwargs={"alpha": 0.6}, ax=ax1)
# plt.plot([0, 60], [0, 60], c=cpal[0], linewidth=3, linestyle='dashed')
plt.plot([0, 15], [0, 15], c=cpal[0], linewidth=3, linestyle='dashed')
plt.grid()
# plt.xlabel("Gaussian Fit $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
# plt.ylabel("FWHM $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
plt.xlabel("CO Integrated Intensity\nfrom Gaussian Fit (K km s$^{-1}$)")
plt.ylabel("CO Integrated Intensity\nin FWHM fit window\n(K km s$^{-1}$)")

ax2 = plt.subplot(122)
hist2d(tab["coldens_HI_gauss"][good_pts] / (inc * hi_mass_conversion).value,
       tab["coldens_HI_FWHM"][good_pts] / (inc * hi_mass_conversion).value,
       data_kwargs={"alpha": 0.6}, ax=ax2)
# plt.plot([0, 25], [0, 25], c=cpal[0], linewidth=3, linestyle='dashed')
plt.plot([0, 2300], [0, 2300], c=cpal[0], linewidth=3, linestyle='dashed')
plt.grid()
# plt.xlabel("Gaussian Fit $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
# plt.ylabel("FWHM $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
plt.xlabel("{\sc HI} Integrated Intensity\nfrom Gaussian Fit (K km s$^{-1}$)")
plt.ylabel("{\sc HI} Integrated Intensity\nin FWHM fit window\n(K km s$^{-1}$)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "coldens_fit_vs_fwhm_check.png"))
plt.savefig(osjoin(fig_path, "coldens_fit_vs_fwhm_check.pdf"))
plt.close()


# Before loading the column densities from the moment arrays, use this
# comparison to find the line width ratios

# Use an error-in-variables approach for finding the ratio. This functions is
# adapted from TurbuStat and fixes the intercept at the origin.

def bayes_linear(x, y, x_err, y_err, nWalkers=10, nBurn=100, nSample=1000,
                 nThin=5, conf_interval=[15.9, 84.1], fix_intercept=False):
    '''
    Fit a line w/ the intercept set to 0.
    '''

    mean_scatter = np.mean(np.sqrt(x_err**2 + y_err**2))
    std_scatter = np.std(np.sqrt(x_err**2 + y_err**2))

    if fix_intercept:
        def _logprob(p, x, y, x_err, y_err):
            theta, var = p[0], p[1]
            if np.abs(theta - np.pi / 4) > np.pi / 4:
                return -np.inf
            if var < 0:
                return -np.inf

            Delta = (np.cos(theta) * y - np.sin(theta) * x)**2
            Sigma = np.sin(theta)**2 * x_err**2 + np.cos(theta)**2 * y_err**2
            lp = -0.5 * np.nansum(Delta / (Sigma + var)) - \
                0.5 * np.nansum(np.log(Sigma + var))
            return lp

        ndim = 2

        p0 = np.zeros((nWalkers, ndim))
        p0[:, 0] = np.tan(np.nanmean(y / x)) + np.random.randn(nWalkers) * 0.1
        p0[:, 1] = np.random.normal(mean_scatter, 0.1 * std_scatter,
                                    size=nWalkers)

    else:
        def _logprob(p, x, y, x_err, y_err):
            theta, bcos, var = p[0], p[1], p[2]
            if np.abs(theta - np.pi / 4) > np.pi / 4:
                return -np.inf
            if var < 0:
                return -np.inf

            Delta = (np.cos(theta) * y - np.sin(theta) * x - bcos)**2
            Sigma = (np.sin(theta))**2 * x_err**2 + \
                (np.cos(theta))**2 * y_err**2
            lp = -0.5 * np.nansum(Delta / (Sigma + var)) - \
                0.5 * np.nansum(np.log(Sigma + var))

            return lp

        ndim = 3
        p0 = np.zeros((nWalkers, ndim))
        p0[:, 0] = np.tan(np.nanmean(y / x)) + np.random.randn(nWalkers) * 0.1
        p0[:, 1] = np.random.randn(nWalkers) * y.std() + y.mean()
        p0[:, 2] = np.random.normal(mean_scatter, 0.1 * std_scatter,
                                    size=nWalkers)

    sampler = emcee.EnsembleSampler(nWalkers, ndim, _logprob,
                                    args=[x, y, x_err, y_err])
    pos, prob, state = sampler.run_mcmc(p0, nBurn)
    sampler.reset()
    sampler.run_mcmc(pos, nSample, thin=nThin)

    if fix_intercept:
        slopes = np.tan(sampler.flatchain[:, 0])
        slope = np.median(slopes)

        var = np.sqrt(sampler.flatchain[:, 1])
        add_stddev = np.median(var)

        params = np.array([slope, add_stddev])
        # Use the percentiles given in conf_interval
        error_intervals = np.empty((2, 2))
        error_intervals[0] = np.percentile(slopes, conf_interval)
        error_intervals[1] = np.percentile(var, conf_interval)

    else:
        slopes = np.tan(sampler.flatchain[:, 0])
        intercepts = sampler.flatchain[:, 1] / np.cos(sampler.flatchain[:, 0])

        slope = np.median(slopes)
        intercept = np.median(intercepts)

        var = np.sqrt(sampler.flatchain[:, 1])
        add_stddev = np.median(var)

        params = np.array([slope, intercept, var])

        # Use the percentiles given in conf_interval
        error_intervals = np.empty((3, 2))
        error_intervals[0] = np.percentile(slopes, conf_interval)
        error_intervals[1] = np.percentile(intercepts, conf_interval)
        error_intervals[2] = np.percentile(var, conf_interval)

    return params, error_intervals, sampler


# Fit an intercept
# params, cis, sampler = \
#     bayes_linear(tab['sigma_HI'][good_pts], tab['sigma_CO'][good_pts],
#                  tab['sigma_stderr_HI'][good_pts],
#                  tab['sigma_stderr_CO'][good_pts],
#                  nBurn=500, nSample=2000, nThin=2)

# Just fit ratio (intercept is 0)
pars_ratio, pars_ratio_ci, sampler_ratio = \
    bayes_linear(tab['sigma_HI'][good_pts], tab['sigma_CO'][good_pts],
                 tab['sigma_stderr_HI'][good_pts],
                 tab['sigma_stderr_CO'][good_pts],
                 nBurn=500, nSample=5000, nThin=1,
                 fix_intercept=True)

slope_ratio = pars_ratio[0]
slope_ratio_ci = pars_ratio_ci[0]

add_stddev_ratio = pars_ratio[1]
add_stddev_ratio_ci = pars_ratio_ci[1]

onecolumn_figure()
hist2d(tab['sigma_HI'][good_pts] / 1000.,
       np.array(tab['sigma_CO'])[good_pts] / 1000., bins=13,
       data_kwargs={"alpha": 0.5})
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")

# slope = params[0]
# inter = params[1] / 1000.
# slope_ci = cis[0]
# inter_cis = cis[1] / 1000.
# plt.plot([4, 12], [4. * slope + inter, 12. * slope + inter],
#          label='Linear Fit')
# plt.fill_between([4, 12], [4. * slope_ci[0] + inter_cis[0],
#                            12. * slope_ci[0] + inter_cis[0]],
#                  [4. * slope_ci[1] + inter_cis[1],
#                   12. * slope_ci[1] + inter_cis[1]],
#                  facecolor=sb.color_palette()[0],
#                  alpha=0.5)
plt.plot([4, 12], [4. * slope_ratio, 12. * slope_ratio],
         '--', color=sb.color_palette()[1], linewidth=3,
         label='Ratio Fit')

# Shade in the region from the added variance
# y_var = add_stddev_ratio * np.sin(np.tan(slope_ratio)) / 1000.

# plt.fill_between([4, 12], [4. * slope_ratio - y_var,
#                            12. * slope_ratio - y_var],
#                  [4. * slope_ratio + y_var,
#                   12. * slope_ratio + y_var],
#                  facecolor=sb.color_palette()[1],
#                  alpha=0.35)

plt.plot([4, 12], [4, 12], '-.', linewidth=3, alpha=0.8,
         color=sb.color_palette()[2],
         label=r'$\sigma_{\rm CO} = \sigma_{\rm HI}$')

plt.axhline(2.6, color=sb.color_palette()[3], linestyle=':',
            alpha=0.5, linewidth=3)

plt.ylim([0.5, 8])
plt.legend(frameon=True, loc='lower right')
plt.grid()

plt.tight_layout()

plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_fit.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_fit.pdf"))
plt.close()

# print("Slope: {0} {1}".format(slope, slope_ci))
# print("Intercept: {0} {1}".format(inter, inter_cis))

print("Ratio Slope: {0} {1}".format(slope_ratio, slope_ratio_ci))
# Ratio Slope: 0.563181915083 [ 0.5621998   0.56414808]

print("Added stddev: {0} {1}".format(add_stddev_ratio, add_stddev_ratio_ci))
# Added stddev: 520.947992994 [ 514.89901222  527.00225878]

# What does this relation look like for line widths from the second moment
co_lwidth = Projection.from_hdu(fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.lwidth.fits"))[0])
hi_lwidth = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['LWidth'])[0])

co_lwidth_vals = co_lwidth.value[tab['ypts'][good_pts], tab['xpts'][good_pts]] / 1000.
hi_lwidth_vals = hi_lwidth.value[tab['ypts'][good_pts], tab['xpts'][good_pts]] / 1000.

# How bad is the relation between the 2nd moment line widths

hist2d(hi_lwidth_vals, co_lwidth_vals, bins=13,
       data_kwargs={"alpha": 0.5})
plt.plot([4, 16], [4. * slope_ratio, 16. * slope_ratio],
         '--', color=sb.color_palette()[1], linewidth=3,
         label='Ratio Fit')

# plt.plot([4, 16], [4. * slope + inter, 16. * slope + inter])
# plt.fill_between([4, 16], [4. * slope_ci[0] + inter_cis[0],
#                            16. * slope_ci[0] + inter_cis[0]],
#                  [4. * slope_ci[1] + inter_cis[1],
#                   16. * slope_ci[1] + inter_cis[1]],
#                  facecolor=sb.color_palette()[0],
#                  alpha=0.5)
plt.plot([4, 16], [4, 16], '--', linewidth=3, alpha=0.8)
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.tight_layout()
plt.grid()
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_2ndmoment.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_2ndmoment.pdf"))
plt.close()

# It's pretty bad

# Peak vs. offset in velocity

hist2d((tab['mean_CO'][good_pts] - tab['mean_HI'][good_pts]) / 1000.,
       tab['amp_CO'][good_pts],
       bins=20,
       data_kwargs={"alpha": 0.5})
plt.grid()
plt.xlabel(r"$|v_{\rm CO} - v_{\rm HI}|$ (km/s)")
plt.ylabel(r"T$_{\rm CO} (K)$")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "vdiff_CO_peak_fit.png"))
plt.savefig(osjoin(fig_path, "vdiff_CO_peak_fit.pdf"))
plt.close()

# Same plotted against HI amplitude

hist2d((tab['mean_CO'][good_pts] - tab['mean_HI'][good_pts]) / 1000.,
       tab['amp_HI'][good_pts],
       bins=20,
       data_kwargs={"alpha": 0.5})
plt.grid()
plt.xlabel(r"$|v_{\rm CO} - v_{\rm HI}|$ (km/s)")
plt.ylabel(r"T$_{\rm HI}$ (K)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "vdiff_HI_peak_fit.png"))
plt.savefig(osjoin(fig_path, "vdiff_HI_peak_fit.pdf"))
plt.close()

# HI vs H2 column density

hist2d(tab['coldens_HI_gauss'][good_pts],
       tab['coldens_CO_gauss'][good_pts],
       bins=30,
       data_kwargs={"alpha": 0.5})
plt.grid()
plt.xlabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")
plt.ylabel(r"$\Sigma_{\rm CO}$ (M$_{\odot}$ pc$^{-2}$)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "Sigma_HI_vs_CO_fit.png"))
plt.savefig(osjoin(fig_path, "Sigma_HI_vs_CO_fit.pdf"))
plt.close()

# Peak temps comparison
hist2d(tab['amp_HI'][good_pts], tab['amp_CO'][good_pts],
       bins=30,
       data_kwargs={"alpha": 0.5})
plt.grid()
plt.xlabel(r"T$_{\rm HI}$ (K)")
plt.ylabel(r"T$_{\rm CO}$ (K)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "peak_HI_vs_CO_fit.png"))
plt.savefig(osjoin(fig_path, "peak_HI_vs_CO_fit.pdf"))
plt.close()

# CO width vs. CO peak
hist2d(tab['sigma_CO'][good_pts] / 1000.,
       tab['amp_CO'][good_pts],
       bins=23,
       data_kwargs={"alpha": 0.5})

# There's an expected correlation if the integrated intensity is constant
onefive_rel = lambda sigma: np.percentile(tab['coldens_CO_gauss'][good_pts] / co21_mass_conversion, 15) / (np.sqrt(2 * np.pi) * sigma)
eightyfive_rel = lambda sigma: np.percentile(tab['coldens_CO_gauss'][good_pts] / co21_mass_conversion, 85) / (np.sqrt(2 * np.pi) * sigma)
median_rel = lambda sigma: np.median(tab['coldens_CO_gauss'][good_pts] / co21_mass_conversion) / (np.sqrt(2 * np.pi) * sigma)

sigma_range = np.linspace(2.5, 8, 100)

plt.plot(sigma_range, median_rel(sigma_range),
         color=sb.color_palette()[0],
         alpha=1.0, linewidth=3,
         linestyle='--')
plt.fill_between(sigma_range, onefive_rel(sigma_range),
                 eightyfive_rel(sigma_range),
                 color=sb.color_palette()[0],
                 alpha=0.3)

plt.grid()
plt.xlabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.ylabel(r"T$_{\rm CO}$ (K)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_CO_vs_peak_CO_fit.png"))
plt.savefig(osjoin(fig_path, "sigma_CO_vs_peak_CO_fit.pdf"))
plt.close()

# HI width vs CO peak

hist2d(tab['sigma_HI'][good_pts] / 1000.,
       tab['amp_CO'][good_pts],
       bins=23,
       data_kwargs={"alpha": 0.5})
plt.grid()
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"T$_{\rm CO}$ (K)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_HI_vs_peak_CO_fit.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_vs_peak_CO_fit.pdf"))
plt.close()

# Load in the column density estimate from the moment arrays
mom_tab = Table.read(fourteenB_HI_data_wGBT_path("tables/column_densities_perpix.fits"))


# Next look at the HI and CO column densities at the same points in the two
# sets. How much does the single component deviate from blindly integrating
# over everything?

ys = tab['ypts'][good_pts]
xs = tab['xpts'][good_pts]

overlap_mask = np.zeros(len(mom_tab["Sigma_Total"]), dtype=bool)

# NOTE!!!: x and y need to be flipped in the moment col dens file!!
for y, x in zip(ys, xs):
    yidx = np.where(mom_tab['xpix'] == y)[0]
    xidx = np.where(mom_tab['ypix'] == x)[0]

    # Look for overlaps
    match = list(set(yidx) & set(xidx))

    if len(match) == 0:
        continue
    if len(match) > 1:
        raise ValueError("Multiple matches? Not possible!")

    overlap_mask[match[0]] = True

twocolumn_twopanel_figure()

fig, axs = plt.subplots(1, 2, sharex=False, sharey=False)

ax1 = axs[1]
# hist2d(mom_tab["Sigma_HI"][overlap_mask],
#        tab['coldens_HI_gauss'][good_pts],
#        data_kwargs={"alpha": 0.6}, ax=ax1)
# Plot everything in data units so we don't need to adopt a X_CO factor for
# paper 2 (not needed for results)
hist2d(mom_tab["Sigma_HI"][overlap_mask] / (inc * hi_mass_conversion).value,
       tab['coldens_HI_gauss'][good_pts] / (inc * hi_mass_conversion).value,
       data_kwargs={"alpha": 0.6}, ax=ax1)
# ax1.plot([0, 25], [0, 25], c=cpal[0], linewidth=3, linestyle='dashed')
ax1.plot([0, 2300], [0, 2300], c=cpal[0], linewidth=3, linestyle='dashed')
ax1.grid()
ax1.set_xlabel("{\sc HI} Integrated Intensity (K km s$^{-1}$)")
ax1.set_ylabel("{\sc HI} Integrated Intensity\n from Gaussian Fit (K km s$^{-1}$)")

# ax1.set_xlabel("Moment $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
# ax1.set_ylabel("Gaussian Fit $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")

ax2 = axs[0]
# hist2d(mom_tab["Sigma_H2"][overlap_mask],
#        tab['coldens_CO_gauss'][good_pts],
#        data_kwargs={"alpha": 0.6}, ax=ax2, bins=20)
hist2d(mom_tab["Sigma_H2"][overlap_mask] / (inc * co21_mass_conversion).value,
       tab['coldens_CO_gauss'][good_pts] / (inc * co21_mass_conversion).value,
       data_kwargs={"alpha": 0.6}, ax=ax2, bins=20)
# ax2.plot([0, 65], [0, 65], c=cpal[0], linewidth=3, linestyle='dashed')
ax2.plot([0, 16], [0, 16], c=cpal[0], linewidth=3, linestyle='dashed')
ax2.grid()
ax2.set_xlabel("CO Integrated Intensity (K km s$^{-1}$)")
ax2.set_ylabel("CO Integrated Intensity\n from Gaussian Fit (K km s$^{-1}$)")

# ax2.set_xlabel("Moment $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
# ax2.set_ylabel("Gaussian Fit $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "coldens_fit_vs_moment_check.png"))
plt.savefig(osjoin(fig_path, "coldens_fit_vs_moment_check.pdf"))
plt.close()

# Does the line width ratio change with galactic radius?
# Compare to the stacked widths
co_radial_fits = Table.read(iram_co21_14B088_data_path("tables/co_hwhm_totalprof_fits_radial_500pc.csv"))
hi_radial_fits = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_radial_500pc.csv"))


rad_bins = np.arange(0, 7.5, 0.5)
med_ratio, bin_edges, cts = \
    binned_statistic(tab['Rgal'][good_pts],
                     tab['sigma_CO'][good_pts] / tab['sigma_HI'][good_pts],
                     bins=rad_bins,
                     statistic=np.median)
bin_cents = (bin_edges[1:] + bin_edges[:-1]) / 2.

ratio_stds, bin_edges, cts = \
    binned_statistic(tab['Rgal'][good_pts],
                     tab['sigma_CO'][good_pts] / tab['sigma_HI'][good_pts],
                     bins=rad_bins,
                     statistic=np.std)

num_in_bins = np.bincount(cts)[1:]
# Num. indep't points divided by number of pixels in one beam.
num_indept = num_in_bins / 41.

onecolumn_figure()

# plt.scatter(tab['Rgal'][good_pts], tab['sigma_CO'][good_pts] / tab['sigma_HI'][good_pts],
#             s=0.1, color='k', rasterized=True, alpha=0.4)

stack_ratio = co_radial_fits['peaksub_sigma'] / hi_radial_fits['peaksub_sigma']
stack_stderr = stack_ratio * \
    np.sqrt((co_radial_fits['peaksub_sigma_up_lim'] /
             co_radial_fits['peaksub_sigma'])**2 +
            (hi_radial_fits['peaksub_sigma_up_lim'] /
             hi_radial_fits['peaksub_sigma'])**2)
plt.errorbar(bin_cents, stack_ratio, yerr=stack_stderr, label='Stack',
             fmt='D-', drawstyle='steps-mid',
             markersize=6)

plt.errorbar(bin_cents, med_ratio, fmt='s-.',
             yerr=ratio_stds / np.sqrt(num_indept), label='Line-of-sight',
             drawstyle='steps-mid', elinewidth=2,
             color=sb.color_palette()[2],
             markersize=6)
# plt.axhline(0.61079832, linestyle='--', color=sb.color_palette()[1])
plt.axhline(slope_ratio, linestyle='--', color=sb.color_palette()[2])

plt.ylabel(r"$\sigma_{\rm CO} / \sigma_{\rm HI}$")
plt.xlabel("Radius (kpc)")

plt.grid()
plt.legend(frameon=True, loc='lower left')

plt.tight_layout()

plt.savefig(osjoin(fig_path, "sigma_ratio_vs_radius.png"))
plt.savefig(osjoin(fig_path, "sigma_ratio_vs_radius.pdf"))
plt.close()

# Now do a fit to the relation per radial bin
# rad_pars = []
# rad_cis = []

# for i, (rad_in, rad_out) in enumerate(zip(rad_bins[:-1], rad_bins[1:])):
#     print(i, rad_in, rad_out)
#     rad_mask = np.logical_and(tab['Rgal'] >= rad_in, tab['Rgal'] < rad_out)
#     rad_par, rad_ci, samp_rad = \
#         bayes_linear(tab['sigma_HI'][good_pts & rad_mask],
#                      tab['sigma_CO'][good_pts & rad_mask],
#                      tab['sigma_stderr_HI'][good_pts & rad_mask],
#                      tab['sigma_stderr_CO'][good_pts & rad_mask],
#                      nBurn=500, nSample=2000, nThin=2)

#     rad_pars.append(rad_par)
#     rad_cis.append(rad_ci)

# rad_pars = np.array(rad_pars)
# rad_cis = np.array(rad_cis)

# plt.errorbar(bin_cents, rad_pars,
#              yerr=[rad_pars - rad_cis[:, 0],
#                    rad_cis[:, 1] - rad_pars],
#              label='Fit', drawstyle='steps-mid')
# plt.errorbar(bin_cents, stack_ratio, yerr=stack_stderr, label='Stack',
#              fmt='D-', drawstyle='steps-mid')

# Is there a pattern in the scatter related to the peak CO temperature?
# Split by quartiles

amp_CO_percs = np.percentile(tab['amp_CO'][good_pts], [0, 25, 50, 75, 100])

twocolumn_figure()

fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)

ratio_ampCO_bins = []

median_CO = []
median_HI = []

samps_ampCO_bins = []

for low, up, ax in zip(amp_CO_percs[:-1], amp_CO_percs[1:], axes.ravel()):

    samples = np.logical_and(tab['amp_CO'] >= low, tab['amp_CO'] <= up)
    samples = np.logical_and(good_pts, samples)

    samps_ampCO_bins.append(samples)

    hist2d(tab['sigma_HI'][samples] / 1000.,
           tab['sigma_CO'][samples] / 1000.,
           bins=10, data_kwargs={"alpha": 0.5},
           ax=ax)
    ax.axvline(np.median(tab['sigma_HI'][samples] / 1000.))
    ax.axhline(np.median(tab['sigma_CO'][samples] / 1000.))

    median_HI.append(np.median(tab['sigma_HI'][samples] / 1000.))
    median_CO.append(np.median(tab['sigma_CO'][samples] / 1000.))

    ax.plot([4, 12], [4, 12], '--', linewidth=3)

    ax.annotate("{0:.2f}--{1:.2f} K".format(low, up),
                xy=(.8, .1), xycoords=ax.transAxes,
                color='k',
                horizontalalignment='center',
                verticalalignment='center',
                bbox={"boxstyle": "square", "facecolor": "w"})

    params_ampCO_bin, cis_ampCO_bin, sampler_ampCO_bin = \
        bayes_linear(tab['sigma_HI'][samples],
                     tab['sigma_CO'][samples],
                     tab['sigma_stderr_HI'][samples],
                     tab['sigma_stderr_CO'][samples],
                     nBurn=500, nSample=1000, nThin=2,
                     fix_intercept=True)
    slope = params_ampCO_bin[0]
    slope_ci = cis_ampCO_bin[0]

    ratio_ampCO_bins.append([slope, slope_ci])

    ax.plot([4, 12], [4. * slope, 12 * slope])
    ax.fill_between([4, 12], [4. * slope_ci[0],
                              12 * slope_ci[0]],
                    [4. * slope_ci[1],
                     12 * slope_ci[1]],
                    facecolor=sb.color_palette()[0],
                    alpha=0.5)

fig.text(0.5, 0.04, r"$\sigma_{\rm HI}$ (km/s)", ha='center')
fig.text(0.04, 0.5, r"$\sigma_{\rm CO}$ (km/s)", va='center',
         rotation='vertical')

fig.savefig(osjoin(fig_path, "sigma_peakCO_bins.png"))
fig.savefig(osjoin(fig_path, "sigma_peakCO_bins.pdf"))
plt.close()

print("Bins: {}".format(amp_CO_percs))
print("Median HI {}".format(median_HI))
print("Median CO {}".format(median_CO))
print("Ratios: {}".format(ratio_ampCO_bins))
# Bins: [ 0.0219169   0.09599498  0.13232244  0.19245095  0.80435075]
# Median HI [7.8417065991420429, 7.5443055678764352, 7.319231497156256, 7.1817473529057017]
# Median CO [5.3499201185891607, 4.4162685854964181, 3.9670922291672444, 3.6825616949564433]
# Ratios: [[0.65052747954370571, array([ 0.64805208,  0.65321125])], [0.58825734786640749, array([ 0.5860826 ,  0.59031582])], [0.54988030892882245, array([ 0.54802902,  0.55170164])], [0.52596686764429645, array([ 0.52458116,  0.52725499])]]

# Where are the different bins located? Different clouds or inter-cloud
# variation?

# hi_mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0])

# twocolumn_figure()
# hi_mom0.quicklook()

# for i, samps in enumerate(samps_ampCO_bins):

#     hi_mom0.FITSFigure._ax1.scatter(tab['xpts'][samps], tab['ypts'][samps],
#                                     c=sb.color_palette()[i])

# You're going to have to stop and zoom in to see the different structure.
# Uncomment to do so.

# They're mostly inter-cloud variation. Most structures have a clear
# progression between the bins from inside out. This also means that
# the line widths in the bins are moderately correlated.
# The larger scatter at lower T_CO bins is probably a S/N

# Compare the HI fit properties w/ and w/o the FWHM mask when fitting
# Fit the same relation as above

params_nomask, cis_nomask, sampler_nomask = \
    bayes_linear(tab['sigma_HI_nomask'][good_pts], tab['sigma_CO'][good_pts],
                 tab['sigma_stderr_HI_nomask'][good_pts],
                 tab['sigma_stderr_CO'][good_pts],
                 nBurn=500, nSample=2000, nThin=2, fix_intercept=True)
ratio = params_nomask[0]
ratio_ci = cis_nomask[0]

print("HI sigma no mask vs. CO")
print("Ratio: {0} {1}".format(ratio, ratio_ci))
# Ratio: 0.504485636162 [ 0.50344384  0.50547394]

onecolumn_figure()

hist2d(tab['sigma_HI'][good_pts] / 1000.,
       tab['sigma_HI_nomask'][good_pts] / 1000.,
       bins=18, data_kwargs={"alpha": 0.5})
plt.plot([4, 12], [4, 12], '--', linewidth=3)

plt.grid()

plt.ylabel(r"No mask $\sigma_{\rm HI}$ (km/s)")
plt.xlabel(r"FWHM mask $\sigma_{\rm HI}$ (km/s)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "sigma_w_wo_masking.png"))
plt.savefig(osjoin(fig_path, "sigma_w_wo_masking.pdf"))
plt.close()

# HI no mask vs CO
hist2d(tab['sigma_HI_nomask'][good_pts] / 1000.,
       tab['sigma_CO'][good_pts] / 1000.,
       bins=18, data_kwargs={"alpha": 0.5})
plt.plot([4, 17.5], [4, 17.5], '-.', linewidth=3, color=sb.color_palette()[2])
plt.plot([4, 17.5], [4. * ratio, 17.5 * ratio],
         color=sb.color_palette()[0])
plt.grid()

plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.xlabel(r"No mask $\sigma_{\rm HI}$ (km/s)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "sigma_HI_nomask_vs_H2_w_fit.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_nomask_vs_H2_w_fit.pdf"))
plt.close()


# Investigate correlations amongst the set of parameters
pan_tab = tab[good_pts].to_pandas()

# Drop some variables
drops = ['RA', 'Dec', 'xpts', 'ypts', 'multicomp_flag_CO', 'multicomp_flag_HI']
drops += [col for col in pan_tab.columns if "nomask" in col]
drops += [col for col in pan_tab.columns if "stderr" in col]
drops += [col for col in pan_tab.columns if "FWHM" in col]
drops += [col for col in pan_tab.columns if "coldens" in col]

drops = list(set(drops))

pan_tab = pan_tab.drop(drops, axis=1)

# Remove _ from names
pan_tab = pan_tab.rename(columns={col: "".join(col.split('_'))
                                  for col in pan_tab.columns})

corr_mat = pan_tab.corr()

mask = np.zeros_like(corr_mat, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True

twocolumn_figure()

he = sb.heatmap(corr_mat, mask=mask)
he.set_xticklabels(he.get_xticklabels(), rotation=90)
he.set_yticklabels(he.get_yticklabels(), rotation=0)

plt.tight_layout()

plt.savefig(osjoin(fig_path, "fit_param_corr_matrix.png"))
plt.savefig(osjoin(fig_path, "fit_param_corr_matrix.pdf"))
plt.close()


def corrfunc(x, y, **kws):
    r, _ = stats.pearsonr(x, y)
    ax = plt.gca()

    if abs(r) >= 0.2:
        ax.annotate("r = {:.2f}".format(r),
                    xy=(.5, .5), xycoords=ax.transAxes,
                    color='r',
                    horizontalalignment='center',
                    verticalalignment='center')
    else:
        ax.annotate("r = {:.2f}".format(r),
                    xy=(.5, .5), xycoords=ax.transAxes,
                    horizontalalignment='center',
                    verticalalignment='center')


def kdeplot_colored(x, y, **kws):
    r, _ = stats.pearsonr(x, y)

    if "cmap" in kws:
        kws.pop('cmap')

    if abs(r) >= 0.2:
        sb.kdeplot(x, y, cmap='Reds_d', **kws)
    else:
        sb.kdeplot(x, y, cmap='Blues_d', **kws)


g = sb.PairGrid(pan_tab, palette=["red"])
g.map_upper(corrfunc)
g.map_diag(sb.distplot, kde=False)
g.map_lower(kdeplot_colored)

plt.tight_layout()

plt.savefig(osjoin(fig_path, "fit_param_pairplot.png"))
plt.savefig(osjoin(fig_path, "fit_param_pairplot.pdf"))
plt.close()

# What fraction of LOS are dominated by H2?
# How does this change from the moments to fits?
h2dom_moments = sum(mom_tab["Ratio"][overlap_mask] >= 1.) / \
    float(overlap_mask.sum())
h2dom_fits = sum((tab['coldens_CO_gauss'] / tab['coldens_HI_gauss'])[good_pts] >= 1.) / \
    float(overlap_mask.sum())
print("Fraction of LOS dominated by H2 from moments: {}".format(h2dom_moments))
# 0.2415
print("Fraction of LOS dominated by H2 from fits: {}".format(h2dom_fits))
# 0.4370

# Median properties?
print("Median Ratio from moments: {}".format(np.median(mom_tab["Ratio"][overlap_mask])))
print("Median Ratio from fits: {}".format(np.median((tab['coldens_CO_gauss'] / tab['coldens_HI_gauss'])[good_pts])))
# Median Ratio from moments: 0.658840902968
# Median Ratio from fits: 0.920772482979

print("Median Gas SD from moments: {}".format(np.median(mom_tab["Sigma_Total"][overlap_mask])))
print("Median Gas SD from fits: {}".format(np.median((tab['coldens_CO_gauss'] +
                                                      tab['coldens_HI_gauss'])[good_pts])))
# Median Gas SD from moments: 20.3851612418
# Median Gas SD from fits: 20.2102707871
