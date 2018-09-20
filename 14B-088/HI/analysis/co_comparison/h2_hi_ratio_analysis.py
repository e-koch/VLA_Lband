

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
from scipy import stats

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, iram_co21_14B088_data_path)
from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_twopanel_figure,
                             twocolumn_figure)

from krumholz_models import krumholz2013_ratio_model, krumholz2013_sigmaHI

cpal = sb.color_palette()

fig_path = allfigs_path("co_vs_hi")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

fig_path_models = allfigs_path("co_vs_hi/h2_formation_models")
if not os.path.exists(fig_path_models):
    os.mkdir(fig_path_models)

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


# Load in the column density estimate from the moment arrays
mom_tab = Table.read(fourteenB_HI_data_wGBT_path("tables/column_densities_perpix.fits"))

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


default_figure()
