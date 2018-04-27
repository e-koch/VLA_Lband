

'''
Compare the H2/HI distribution to theory.

The outputs from co_hi_linewidth_ratio.py and h2_hi_ratios.py should be
available. The former finds the column densities with a single Gaussian fit,
and the latter uses the moment arrays.

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

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, iram_co21_14B088_data_path)
from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_twopanel_figure)

from krumholz_models import krumholz2013_ratio_model, krumholz2013_sigmaHI

fig_path = allfigs_path("co_vs_hi")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

fig_path_models = allfigs_path("co_vs_hi/h2_formation_models")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

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
twocolumn_twopanel_figure()
cpal = sb.color_palette()

ax1 = plt.subplot(121)

hist2d(tab["coldens_CO_FWHM"][good_pts],
       tab["coldens_CO_gauss"][good_pts],
       data_kwargs={"alpha": 0.6}, ax=ax1)
plt.plot([0, 60], [0, 60], c=cpal[0], linewidth=3, linestyle='dashed')
plt.grid()
plt.xlabel("Gaussian Fit $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
plt.ylabel("FWHM $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")

ax2 = plt.subplot(122)
hist2d(tab["coldens_HI_FWHM"][good_pts],
       tab["coldens_HI_gauss"][good_pts],
       data_kwargs={"alpha": 0.6}, ax=ax2)
plt.plot([0, 25], [0, 25], c=cpal[0], linewidth=3, linestyle='dashed')
plt.grid()
plt.xlabel("Gaussian Fit $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
plt.ylabel("FWHM $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")

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

    if fix_intercept:
        def _logprob(p, x, y, x_err, y_err):
            theta = p[0]
            if np.abs(theta - np.pi / 4) > np.pi / 4:
                return -np.inf

            Delta = (np.cos(theta) * y - np.sin(theta) * x)**2
            Sigma = np.sin(theta)**2 * x_err**2 + np.cos(theta)**2 * y_err**2
            lp = -0.5 * np.nansum(Delta / Sigma) - \
                0.5 * np.nansum(np.log(Sigma))
            return lp

        ndim = 2

        p0 = np.zeros((nWalkers, ndim))
        p0[:, 0] = np.pi / 4 + np.random.randn(nWalkers) * 0.1

    else:
        def _logprob(p, x, y, x_err, y_err):
            theta, b = p[0], p[1]
            if np.abs(theta - np.pi / 4) > np.pi / 4:
                return -np.inf
            Delta = (np.cos(theta) * y - np.sin(theta) * x - b * np.cos(theta))**2
            Sigma = (np.sin(theta))**2 * x_err**2 + (np.cos(theta))**2 * y_err**2
            lp = -0.5 * np.nansum(Delta / Sigma) - 0.5 * np.nansum(np.log(Sigma))

            return lp

        ndim = 2
        p0 = np.zeros((nWalkers, ndim))
        p0[:, 0] = np.pi / 4 + np.random.randn(nWalkers) * 0.1
        p0[:, 1] = np.random.randn(nWalkers) * y.std() + y.mean()

    sampler = emcee.EnsembleSampler(nWalkers, ndim, _logprob,
                                    args=[x, y, x_err, y_err])
    pos, prob, state = sampler.run_mcmc(p0, nBurn)
    sampler.reset()
    sampler.run_mcmc(pos, nSample, thin=nThin)

    if fix_intercept:
        slopes = np.tan(sampler.flatchain[:, 0])
        slope = np.median(slopes)
        params = np.array([slope])
        # Use the percentiles given in conf_interval
        error_intervals = np.percentile(slopes, conf_interval)

    else:
        slopes = np.tan(sampler.flatchain[:, 0])
        intercepts = sampler.flatchain[:, 1]

        slope = np.median(slopes)
        intercept = np.median(intercepts)

        params = np.array([slope, intercept])

        # Use the percentiles given in conf_interval
        error_intervals = np.empty((2, 2))
        error_intervals[0] = np.percentile(slopes, conf_interval)
        error_intervals[1] = np.percentile(intercepts, conf_interval)

    return params, error_intervals, sampler


# Fit an intercept
params, cis, sampler = \
    bayes_linear(tab['sigma_HI'][good_pts], tab['sigma_CO'][good_pts],
                 tab['sigma_stderr_HI'][good_pts],
                 tab['sigma_stderr_CO'][good_pts],
                 nBurn=500, nSample=2000, nThin=2)

# Just fit ratio (intercept is 0)
slope_ratio, slope_ratio_ci, sampler_ratio = \
    bayes_linear(tab['sigma_HI'][good_pts], tab['sigma_CO'][good_pts],
                 tab['sigma_stderr_HI'][good_pts],
                 tab['sigma_stderr_CO'][good_pts],
                 nBurn=500, nSample=2000, nThin=2,
                 fix_intercept=True)

onecolumn_figure()
hist2d(tab['sigma_HI'][good_pts] / 1000.,
       np.array(tab['sigma_CO'])[good_pts] / 1000., bins=13,
       data_kwargs={"alpha": 0.5})
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")

slope = params[0]
inter = params[1] / 1000.
slope_ci = cis[0]
inter_cis = cis[1] / 1000.
plt.plot([4, 12], [4. * slope + inter, 12. * slope + inter],
         label='Linear Fit')
plt.fill_between([4, 12], [4. * slope_ci[0] + inter_cis[0],
                           12. * slope_ci[0] + inter_cis[0]],
                 [4. * slope_ci[1] + inter_cis[1],
                  12. * slope_ci[1] + inter_cis[1]],
                 facecolor=sb.color_palette()[0],
                 alpha=0.5)
plt.plot([4, 12], [4. * slope_ratio[0], 12. * slope_ratio[0]],
         '--', color=sb.color_palette()[1], linewidth=3,
         label='Ratio Fit')

plt.plot([4, 12], [4, 12], '-.', linewidth=3, alpha=0.8,
         color=sb.color_palette()[2],
         label=r'$\sigma_{\rm CO} = \sigma_{\rm HI}$')

plt.axhline(2.6, color=sb.color_palette()[3], linestyle=':',
            alpha=0.5, linewidth=3)

plt.ylim([0.5, 8])
plt.legend(frameon=True)

plt.tight_layout()

plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_fit.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_fit.pdf"))
plt.close()

print("Slope: {0} {1}".format(slope, slope_ci))
print("Intercept: {0} {1}".format(inter, inter_cis))
# Slope: 0.845309187275 [ 0.83697267  0.85358567]
# Intercept: -1.82702574095 [-1.89175737 -1.76184693]
print("Ratio Slope: {0} {1}".format(slope_ratio, slope_ratio_ci))
# Ratio Slope: [ 0.6108361] [ 0.60905739  0.61261014]

# What does this relation look like for line widths from the second moment
co_lwidth = Projection.from_hdu(fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.lwidth.fits"))[0])
hi_lwidth = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['LWidth'])[0])

co_lwidth_vals = co_lwidth.value[tab['ypts'][good_pts], tab['xpts'][good_pts]] / 1000.
hi_lwidth_vals = hi_lwidth.value[tab['ypts'][good_pts], tab['xpts'][good_pts]] / 1000.

# How bad is the relation between the 2nd moment line widths

hist2d(hi_lwidth_vals, co_lwidth_vals, bins=13,
       data_kwargs={"alpha": 0.5})
plt.plot([4, 16], [4. * slope + inter, 16. * slope + inter])
plt.fill_between([4, 16], [4. * slope_ci[0] + inter_cis[0],
                           16. * slope_ci[0] + inter_cis[0]],
                 [4. * slope_ci[1] + inter_cis[1],
                  16. * slope_ci[1] + inter_cis[1]],
                 facecolor=sb.color_palette()[0],
                 alpha=0.5)
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

hist2d(tab['coldens_HI_gauss'][good_pts], tab['coldens_CO_gauss'][good_pts],
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
hist2d(mom_tab["Sigma_HI"][overlap_mask],
       tab['coldens_HI_gauss'][good_pts],
       data_kwargs={"alpha": 0.6}, ax=ax1)
ax1.plot([0, 25], [0, 25], c=cpal[0], linewidth=3, linestyle='dashed')
ax1.grid()
ax1.set_xlabel("Moment $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")
ax1.set_ylabel("Gaussian Fit $\Sigma_{\mathrm{HI}}$ (M$_{\odot}$ pc$^{-2}$)")

ax2 = axs[0]
hist2d(mom_tab["Sigma_H2"][overlap_mask],
       tab['coldens_CO_gauss'][good_pts],
       data_kwargs={"alpha": 0.6}, ax=ax2, bins=20)
ax2.plot([0, 65], [0, 65], c=cpal[0], linewidth=3, linestyle='dashed')
ax2.grid()
ax2.set_xlabel("Moment $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
ax2.set_ylabel("Gaussian Fit $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "coldens_fit_vs_moment_check.png"))
plt.savefig(osjoin(fig_path, "coldens_fit_vs_moment_check.pdf"))
plt.close()

# Create the log ratio vs. total Sigma plots
twocolumn_twopanel_figure()

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)

ax1 = axs[1]
hist2d((tab['coldens_HI_gauss'] + tab['coldens_CO_gauss'])[good_pts],
       np.log10((tab['coldens_CO_gauss'] / tab['coldens_HI_gauss'])[good_pts]),
       bins=18, data_kwargs={"alpha": 0.5},
       ax=ax1)
ax1.set_xlabel(r"$\Sigma_{\mathrm{Total}}$ (M$_{\odot}$ pc$^{-2}$)")

# Overplot the Krumholz 2013 model
sigma_t = np.linspace(5, 75, 100)
ax1.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=1, Z=0.5)), "-",
         label="c=1, Z=0.5")
ax1.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=1, Z=1.0)), "--",
         label="c=1, Z=1.0")
ax1.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=2, Z=0.5)), "-.",
         label="c=2, Z=0.5")
ax1.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=3, Z=0.5)), ":",
         label="c=3, Z=0.5")
ax1.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=3, Z=1.0)), "-",
         label="c=3, Z=1.0")
ax1.set_ylim([-1.2, 0.9])
ax1.grid()

ax1.annotate("Gaussian Fit", (0.15, 0.88),
             xycoords='axes fraction', color='k',
             fontsize=12,
             bbox={"boxstyle": "square", "facecolor": "w"})

ax2 = axs[0]
hist2d(mom_tab["Sigma_Total"][overlap_mask],
       np.log10(mom_tab["Ratio"][overlap_mask]),
       data_kwargs={"alpha": 0.6},
       ax=ax2, label='_nolegend_')
ax2.set_xlabel(r"$\Sigma_{\mathrm{Total}}$ (M$_{\odot}$ pc$^{-2}$)")
ax2.set_ylabel(r"log $\Sigma_{\mathrm{H2}} / \Sigma_{\mathrm{HI}}$")

ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=1, Z=0.5)), "-",
         label="c=1, Z=0.5")
ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=1, Z=1.0)), "--",
         label="c=1, Z=1.0")
ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=2, Z=0.5)), "-.",
         label="c=2, Z=0.5")
ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=3, Z=0.5)), ":",
         label="c=3, Z=0.5")
ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=3, Z=1.0)), "-",
         label="c=3, Z=1.0")

ax2.grid()
ax2.set_ylim([-1.2, 1.0])
ax2.set_xlim([5, 80])

ax2.annotate("Moment", (0.15, 0.88),
             xycoords='axes fraction', color='k',
             fontsize=12,
             bbox={"boxstyle": "square", "facecolor": "w"})
# For some reason the 2D histogram is picking up on a column name...
handles, labels = plt.gca().get_legend_handles_labels()
ax2.legend(handles[1:], labels[1:], loc='lower right', frameon=True)

plt.tight_layout()
save_name = "ratio_totalsigma_w_krumholzmodel_perpix_feather_moment_vs_fit"
plt.savefig(osjoin(fig_path_models, "{}.pdf".format(save_name)))
plt.savefig(osjoin(fig_path_models, "{}.png".format(save_name)))
plt.close()


# Sigma HI vs Sigma Total

twocolumn_twopanel_figure()

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)

ax1 = axs[1]
hist2d((tab['coldens_HI_gauss'] + tab['coldens_CO_gauss'])[good_pts],
       tab['coldens_HI_gauss'][good_pts],
       bins=18, data_kwargs={"alpha": 0.5},
       ax=ax1)
ax1.set_xlabel(r"$\Sigma_{\mathrm{Total}}$ (M$_{\odot}$ pc$^{-2}$)")

# Overplot the Krumholz 2013 model
ax1.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=1, Z=1.0), "-.",
         label="c=1, Z=1.0", linewidth=2, alpha=0.95,
         color=cpal[0])
ax1.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=2, Z=0.5), ":",
         label="c=2, Z=0.5", linewidth=2, alpha=0.95,
         color=cpal[1])
ax1.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=3, Z=0.5), "-",
         label="c=3, Z=0.5", linewidth=2, alpha=0.95,
         color=cpal[2])
ax1.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=3, Z=1.0), "--",
         label="c=3, Z=1.0", linewidth=2, alpha=0.95,
         color=cpal[3])

ax1.plot([5, 26], [5, 26], '-', linewidth=4, alpha=0.6, color=cpal[5])

ax1.grid()

ax1.annotate("Gaussian Fit", (0.7, 0.88),
             xycoords='axes fraction', color='k',
             fontsize=12,
             bbox={"boxstyle": "square", "facecolor": "w"})

ax2 = axs[0]
hist2d(mom_tab["Sigma_Total"][overlap_mask],
       mom_tab["Sigma_HI"][overlap_mask],
       data_kwargs={"alpha": 0.6},
       ax=ax2, label='_nolegend_')
ax2.set_xlabel(r"$\Sigma_{\mathrm{Total}}$ (M$_{\odot}$ pc$^{-2}$)")
ax2.set_ylabel(r"$\Sigma_{\mathrm{HI}}$")

ax2.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=1, Z=1.0), "-.",
         label="c=1, Z=1.0", linewidth=2, alpha=0.95,
         color=cpal[0])
ax2.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=2, Z=0.5), ":",
         label="c=2, Z=0.5", linewidth=2, alpha=0.95,
         color=cpal[1])
ax2.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=3, Z=0.5), "-",
         label="c=3, Z=0.5", linewidth=2, alpha=0.95,
         color=cpal[2])
ax2.plot(sigma_t, krumholz2013_sigmaHI(sigma_t, c=3, Z=1.0), "--",
         label="c=3, Z=1.0", linewidth=2, alpha=0.95,
         color=cpal[3])


ax2.plot([5, 26], [5, 26], '-', linewidth=4, alpha=0.6, color=cpal[5])

ax2.grid()
ax2.set_ylim([-2, 26])
ax2.set_xlim([5, 75])

ax2.annotate("Moment", (0.79, 0.88),
             xycoords='axes fraction', color='k',
             fontsize=12,
             bbox={"boxstyle": "square", "facecolor": "w"})
# For some reason the 2D histogram is picking up on a column name...
handles, labels = plt.gca().get_legend_handles_labels()
ax2.legend(handles[1:], labels[1:], loc='lower right', frameon=True, ncol=2)

plt.tight_layout()

save_name = "sigmahi_totalsigma_w_krumholzmodel_perpix_feather_moment_vs_fit"
plt.savefig(osjoin(fig_path_models, "{}.pdf".format(save_name)))
plt.savefig(osjoin(fig_path_models, "{}.png".format(save_name)))
plt.close()


# What fraction of LOS are dominated by H2?
# How does this change from the moments to fits?
h2dom_moments = sum(mom_tab["Ratio"][overlap_mask] >= 1.) / \
    float(overlap_mask.sum())
h2dom_fits = sum((tab['coldens_CO_gauss'] / tab['coldens_HI_gauss'])[good_pts] >= 1.) / \
    float(overlap_mask.sum())
print("Fraction of LOS dominated by H2 from moments: {}".format(h2dom_moments))
# 0.2369
print("Fraction of LOS dominated by H2 from fits: {}".format(h2dom_fits))
# 0.43265

# Median properties?
print("Median Ratio from moments: {}".format(np.median(mom_tab["Ratio"][overlap_mask])))
print("Median Ratio from fits: {}".format(np.median((tab['coldens_CO_gauss'] / tab['coldens_HI_gauss'])[good_pts])))
# Median Ratio from moments: 0.661082448156
# Median Ratio from fits: 0.917707557645

print("Median Gas SD from moments: {}".format(np.median(mom_tab["Sigma_Total"][overlap_mask])))
print("Median Gas SD from fits: {}".format(np.median((tab['coldens_CO_gauss'] + tab['coldens_HI_gauss'])[good_pts])))
# Median Gas SD from moments: 20.3874723999
# Median Gas SD from fits: 20.1290850074

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

lower_ratio, bin_edges, cts = \
    binned_statistic(tab['Rgal'][good_pts],
                     tab['sigma_CO'][good_pts] / tab['sigma_HI'][good_pts],
                     bins=rad_bins,
                     statistic=lambda x: np.percentile(x, 15))

upper_ratio, bin_edges, cts = \
    binned_statistic(tab['Rgal'][good_pts],
                     tab['sigma_CO'][good_pts] / tab['sigma_HI'][good_pts],
                     bins=rad_bins,
                     statistic=lambda x: np.percentile(x, 85))

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
             fmt='D-', drawstyle='steps-mid')

plt.errorbar(bin_cents, med_ratio, fmt='o-.',
             yerr=[med_ratio - lower_ratio,
                   upper_ratio - med_ratio], label='Fit',
             drawstyle='steps-mid')
# plt.axhline(0.61079832, linestyle='--', color=sb.color_palette()[2])
plt.axhline(slope_ratio, linestyle='--', color=sb.color_palette()[2])

plt.ylabel(r"$\sigma_{\rm CO} / \sigma_{\rm HI}$")
plt.xlabel("Radius (kpc)")

plt.grid()
plt.legend(frameon=True, loc='upper right')

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

# Compare the HI fit properties w/ and w/o the FWHM mask when fitting
# Fit the same relation as above

params_nomask, cis_nomask, sampler_nomask = \
    bayes_linear(tab['sigma_HI_nomask'][good_pts], tab['sigma_CO'][good_pts],
                 tab['sigma_stderr_HI_nomask'][good_pts],
                 tab['sigma_stderr_CO'][good_pts],
                 nBurn=500, nSample=2000, nThin=2)
slope = params_nomask[0]
inter = params_nomask[1] / 1000.
slope_ci = cis_nomask[0]
inter_cis = cis_nomask[1] / 1000.

print("HI sigma no mask vs. CO")
print("Slope: {0} {1}".format(slope, slope_ci))
print("Intercept: {0} {1}".format(inter, inter_cis))
# Slope: 0.521911603815 [ 0.51385501  0.52925186]
# Intercept: 0.149811659879 [ 0.08597853  0.21841618]

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
plt.plot([4, 17.5], [4, 17.5], '--', linewidth=3)
plt.plot([4, 17.5], [4. * slope + inter, 17.5 * slope + inter])
plt.fill_between([4, 17.5], [4. * slope_ci[0] + inter_cis[0],
                             17.5 * slope_ci[0] + inter_cis[0]],
                 [4. * slope_ci[1] + inter_cis[1],
                  17.5 * slope_ci[1] + inter_cis[1]],
                 facecolor=sb.color_palette()[0],
                 alpha=0.5)

plt.grid()

plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.xlabel(r"No mask $\sigma_{\rm HI}$ (km/s)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "sigma_HI_nomask_vs_H2_w_fit.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_nomask_vs_H2_w_fit.pdf"))
plt.close()


default_figure()
