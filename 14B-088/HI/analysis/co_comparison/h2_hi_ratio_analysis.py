

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

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, iram_co21_14B088_data_path)
from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_twopanel_figure)

from krumholz_models import krumholz2013_ratio_model, krumholz_maxhi_sigma

fig_path = allfigs_path("co_vs_hi")
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

# Show that the Gaussian model is a decent representation
twocolumn_twopanel_figure()
cpal = sb.color_palette()

ax1 = plt.subplot(121)

hist2d(tab["coldens_CO_FWHM"][good_pts],
       tab["coldens_CO_gauss"][good_pts],
       data_kwargs={"alpha": 0.6}, ax=ax1)
plt.plot([0, 50], [0, 50], c=cpal[0], linewidth=3, linestyle='dashed')
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
                 nThin=5,
                 conf_interval=[15.9, 84.1]):
    '''
    Fit a line w/ the intercept set to 0.
    '''

    # def _logprob(p, x, y, x_err, y_err):
    #     theta = p[0]
    #     if np.abs(theta - np.pi / 4) > np.pi / 4:
    #         return -np.inf

    #     Delta = (np.cos(theta) * y - np.sin(theta) * x)**2
    #     Sigma = np.sin(theta)**2 * x_err**2 + np.cos(theta)**2 * y_err**2
    #     lp = -0.5 * np.nansum(Delta / Sigma) - 0.5 * np.nansum(np.log(Sigma))
    #     return lp

    # ndim = 2

    # p0 = np.zeros((nWalkers, ndim))
    # p0[:, 0] = np.pi / 4 + np.random.randn(nWalkers) * 0.1

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
    # slopes = np.tan(sampler.flatchain[:, 0])
    # slope = np.median(slopes)
    # # Use the percentiles given in conf_interval
    # error_interval = np.percentile(slopes, conf_interval)

    slopes = np.tan(sampler.flatchain[:, 0])
    intercepts = sampler.flatchain[:, 1]

    slope = np.median(slopes)
    intercept = np.median(intercepts)

    params = np.array([slope, intercept])

    # Use the percentiles given in conf_interval
    error_intervals = np.empty((2, 2))
    error_intervals[0] = np.percentile(slopes, conf_interval)
    error_intervals[1] = np.percentile(intercepts, conf_interval)

    # return slope, error_interval, sampler
    return params, error_intervals, sampler


# slope, slope_ci, sampler = \
params, cis, sampler = \
    bayes_linear(tab['sigma_HI'][good_pts], tab['sigma_CO'][good_pts],
                 tab['sigma_stderr_HI'][good_pts],
                 tab['sigma_stderr_CO'][good_pts],
                 nBurn=1000, nSample=5000, nThin=2)

onecolumn_figure()
hist2d(tab['sigma_HI'][good_pts] / 1000.,
       tab['sigma_CO'][good_pts] / 1000., bins=13,
       data_kwargs={"alpha": 0.5})
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")

slope = params[0]
inter = params[1] / 1000.
slope_ci = cis[0]
inter_cis = cis[1] / 1000.
plt.plot([4, 12], [4. * slope + inter, 12. * slope + inter])
plt.fill_between([4, 12], [4. * slope_ci[0] + inter_cis[0],
                           12. * slope_ci[0] + inter_cis[0]],
                 [4. * slope_ci[1] + inter_cis[1],
                  12. * slope_ci[1] + inter_cis[1]],
                 facecolor=sb.color_palette()[0],
                 alpha=0.5)
plt.plot([4, 12], [4, 12], '--', linewidth=3, alpha=0.8)
# plt.fill_between([4, 12], [4. * slope_ci[0], 12. * slope_ci[0]],
#                  [4. * slope_ci[1], 12. * slope_ci[1]], facecolor='gray',
#                  alpha=0.5)

plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_fit.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_fit.pdf"))
plt.close()

print("Slope: {0} {1}".format(slope, slope_ci))
print("Intercept: {0} {1}".format(inter, inter_cis))
# Slope: 0.845577294464 [ 0.83733223  0.8538308 ]
# Intercept: -1.82822168673 [-1.89184967 -1.76575023]

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
ax2.plot([0, 50], [0, 50], c=cpal[0], linewidth=3, linestyle='dashed')
ax2.grid()
ax2.set_xlabel("Moment $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
ax2.set_ylabel("Gaussian Fit $\Sigma_{\mathrm{H2}}$ (M$_{\odot}$ pc$^{-2}$)")
ax2.set_ylim([3, 50])
ax2.set_xlim([3, 50])

plt.tight_layout()

plt.savefig(osjoin(fig_path, "coldens_fit_vs_moment_check.png"))
plt.savefig(osjoin(fig_path, "coldens_fit_vs_moment_check.pdf"))
plt.close()

print(argh)

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
sigma_t = np.linspace(5, 70, 100)
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

ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=1, Z=1.0)), "--",
         label="c=1, Z=1.0")
ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=2, Z=0.5)), "-.",
         label="c=2, Z=0.5")
ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=3, Z=0.5)), ":",
         label="c=3, Z=0.5")
ax2.plot(sigma_t, np.log10(krumholz2013_ratio_model(sigma_t, c=3, Z=1.0)), "-",
         label="c=3, Z=1.0")

ax2.grid()
ax2.set_ylim([-1.2, 0.9])

ax2.annotate("Moment", (0.15, 0.88),
             xycoords='axes fraction', color='k',
             fontsize=12,
             bbox={"boxstyle": "square", "facecolor": "w"})
# For some reason the 2D histogram is picking up on a column name...
handles, labels = plt.gca().get_legend_handles_labels()
ax2.legend(handles[1:], labels[1:], loc='lower right', frameon=True)

plt.tight_layout()

save_name = "ratio_totalsigma_w_krumholzmodel_perpix_feather_moment_vs_fit"
plt.savefig(osjoin(fig_path, "{}.pdf".format(save_name)))
plt.savefig(osjoin(fig_path, "{}.png".format(save_name)))
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
ax1.axhline(krumholz_maxhi_sigma(Z=1.0, c=1).value,
            linestyle='-.', label='Z=1.0, c=1',
            c=cpal[0])
ax1.axhline(krumholz_maxhi_sigma(Z=0.5, c=2).value,
            linestyle=':', label='Z=0.5, c=2',
            c=cpal[1])
ax1.axhline(krumholz_maxhi_sigma(Z=0.5, c=3).value,
            linestyle='-', label='Z=0.5, c=3',
            c=cpal[2])
ax1.axhline(krumholz_maxhi_sigma(Z=1., c=3.).value,
            linestyle='--', label='Z=1.0, c=3',
            c=cpal[3])

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
ax2.set_ylabel(r"log $\Sigma_{\mathrm{H2}} / \Sigma_{\mathrm{HI}}$")


ax2.axhline(krumholz_maxhi_sigma(Z=1.0, c=1).value,
            linestyle='-.', label='Z=1.0, c=1',
            c=cpal[0])
ax2.axhline(krumholz_maxhi_sigma(Z=0.5, c=2).value,
            linestyle=':', label='Z=0.5, c=2',
            c=cpal[1])
ax2.axhline(krumholz_maxhi_sigma(Z=0.5, c=3).value,
            linestyle='-', label='Z=0.5, c=3',
            c=cpal[2])
ax2.axhline(krumholz_maxhi_sigma(Z=1., c=3.).value,
            linestyle='--', label='Z=1.0, c=3',
            c=cpal[3])

ax2.grid()
ax2.set_ylim([-2, 26])

ax2.annotate("Moment", (0.79, 0.88),
             xycoords='axes fraction', color='k',
             fontsize=12,
             bbox={"boxstyle": "square", "facecolor": "w"})
# For some reason the 2D histogram is picking up on a column name...
handles, labels = plt.gca().get_legend_handles_labels()
ax2.legend(handles[1:], labels[1:], loc='lower right', frameon=True, ncol=2)

plt.tight_layout()

save_name = "sigmahi_totalsigma_w_krumholzmodel_perpix_feather_moment_vs_fit"
plt.savefig(osjoin(fig_path, "{}.pdf".format(save_name)))
plt.savefig(osjoin(fig_path, "{}.png".format(save_name)))
plt.close()

default_figure()
