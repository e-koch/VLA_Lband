
'''
Regrid the Corbelli+2017 cloud catalogue asgn file
to the 14B-088
'''

import numpy as np
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import os
from os.path import join as osjoin
from corner import hist2d, corner
import emcee
import matplotlib.pyplot as plt
import seaborn as sb

from paths import (fourteenB_HI_data_wGBT_path,
                   data_path, allfigs_path,
                   iram_co21_14B088_data_path)
from plotting_styles import (default_figure,
                             twocolumn_figure,
                             twocolumn_twopanel_figure,
                             onecolumn_figure)
from galaxy_params import gal_feath as gal


fig_path = allfigs_path("co_vs_hi")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

# Load GMC catalogue from Corbelli+17

gmc_tab = Table.read(osjoin(data_path,
                            'Corbelli_17_catalogues',
                            'J_A+A_601_A146_table5.dat.fits'))

# Load in the Gaussian HI and CO fit table.
tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_co_gaussfit_column_densities_perpix.fits"))

# Don't consider the "bad fits" that are probably due to multiple components
good_pts = np.logical_and(~tab['multicomp_flag_HI'],
                          ~tab['multicomp_flag_CO'])
good_pts = np.logical_and(good_pts,
                          tab["sigma_HI"] > 3800)
# Minimum CO line width of one channel.
good_pts = np.logical_and(good_pts,
                          tab["sigma_CO"] >= 2600)

# Assign the cloud type based on the closest GMC location in the catalogue

cloud_posns = SkyCoord(gmc_tab['RAdeg'], gmc_tab['DEdeg'])
pix_posns = SkyCoord(tab['RA'], tab['Dec'])

# Just do this for every pixel, even if it's not in the good mask
cloud_types = []
for i in range(len(tab)):

    cloud_idx = pix_posns[i].separation(cloud_posns).argmin()

    cloud_types.append(gmc_tab['Type'][cloud_idx])

cloud_types = Column(cloud_types)

tab.add_column(cloud_types, name='cloud_type')


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


twocolumn_figure()

# Original included the "D" column. But since it's poorly
# sampled, remove from plotting.

# fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

fitted_ratios = []
fitted_ratio_errs = []

intrinsic_scatter = []
intrinsic_scatter_errs = []

# for c_type, ax in zip(['A', 'B', 'C', 'D'], axs.ravel()):
for c_type, ax in zip(['A', 'B', 'C'], axs.ravel()):

    type_pts = np.logical_and(good_pts,
                              tab['cloud_type'] == c_type)

    pars_ratio, pars_ratio_ci, sampler_ratio = \
        bayes_linear(tab['sigma_HI'][type_pts], tab['sigma_CO'][type_pts],
                     tab['sigma_stderr_HI'][type_pts],
                     tab['sigma_stderr_CO'][type_pts],
                     nBurn=500, nSample=2000, nThin=1,
                     fix_intercept=True)

    slope_ratio = pars_ratio[0]
    slope_ratio_ci = pars_ratio_ci[0]

    fitted_ratios.append(slope_ratio)
    fitted_ratio_errs.append(slope_ratio_ci)

    add_stddev_ratio = pars_ratio[1]
    add_stddev_ratio_ci = pars_ratio_ci[1]

    intrinsic_scatter.append(add_stddev_ratio)
    intrinsic_scatter_errs.append(add_stddev_ratio_ci)

    hist2d(tab['sigma_HI'][type_pts] / 1000.,
           np.array(tab['sigma_CO'])[type_pts] / 1000.,
           bins=10,
           data_kwargs={"alpha": 0.5},
           ax=ax, range=[(2, 12), (2, 12)])

    ax.plot([2, 12], [2. * slope_ratio, 12. * slope_ratio],
            '--', color=sb.color_palette()[2], linewidth=3,
            label='Ratio Fit')

    ax.plot([2, 12], [2, 12], '-', linewidth=3, alpha=0.8,
            color=sb.color_palette()[3],
            label=r'$\sigma_{\rm CO} = \sigma_{\rm HI}$')

    ax.text(0.1, 0.8, c_type, transform=ax.transAxes, fontsize=20,
            bbox={"boxstyle": "round", "facecolor": "w"})

    ax.axhline(2.6, color=sb.color_palette()[4], linestyle=':',
               alpha=0.5, linewidth=3)
    ax.axhline(8.0, color=sb.color_palette()[4], linestyle=':',
               alpha=0.5, linewidth=3)

axs[1].set_xlabel(r"$\sigma_{\rm HI}$ (km/s)")
axs[0].set_ylabel(r"$\sigma_{\rm CO}$ (km/s)")

# axs[1, 0].set_xlabel(r"$\sigma_{\rm HI}$ (km/s)")
# axs[1, 1].set_xlabel(r"$\sigma_{\rm HI}$ (km/s)")
# axs[0, 0].set_ylabel(r"$\sigma_{\rm CO}$ (km/s)")
# axs[1, 0].set_ylabel(r"$\sigma_{\rm CO}$ (km/s)")

plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_cloudtype.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_cloudtype.pdf"))
plt.close()

# for i, c_type in enumerate(['A', 'B', 'C', 'D']):
for i, c_type in enumerate(['A', 'B', 'C']):
    print("Ratio of {0}: {1}  {2}".format(c_type, fitted_ratios[i],
                                          fitted_ratio_errs[i]))
    print("Scatter of {0}: {1}  {2}".format(c_type, intrinsic_scatter[i],
                                            intrinsic_scatter_errs[i]))

# Ratio of A: 0.6097540942545268  [0.60628682 0.61317084]
# Scatter of A: 586.7949472028565  [566.34715418 607.77593667]
# Ratio of B: 0.5976651855038639  [0.59490276 0.60043233]
# Scatter of B: 597.2436450732399  [580.43501309 614.19493902]
# Ratio of C: 0.5500753222800747  [0.54896911 0.55114205]
# Scatter of C: 486.39031452483346  [479.67843996 492.94783184]
# Ratio of D: 0.5237602437236004  [0.52091347 0.52662596]
# Scatter of D: 158.28200242599462  [129.2197323  188.37514809]

# Make a version w/o the fits

twocolumn_twopanel_figure()

fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

for c_type, ax in zip(['A', 'B', 'C'], axs.ravel()):

    type_pts = np.logical_and(good_pts,
                              tab['cloud_type'] == c_type)

    hist2d(tab['sigma_HI'][type_pts] / 1000.,
           np.array(tab['sigma_CO'])[type_pts] / 1000., bins=10,
           data_kwargs={"alpha": 0.5},
           ax=ax, range=[(2, 12), (2, 12)])

    ax.plot([2, 12], [2, 12], '-', linewidth=3, alpha=0.8,
            color=sb.color_palette()[3],
            label=r'$\sigma_{\rm CO} = \sigma_{\rm HI}$')

    ax.text(0.1, 0.8, c_type, transform=ax.transAxes, fontsize=20,
            bbox={"boxstyle": "round", "facecolor": "w"})

    ax.axhline(2.6, color=sb.color_palette()[4], linestyle=':',
               alpha=0.5, linewidth=3)
    ax.axhline(8.0, color=sb.color_palette()[4], linestyle=':',
               alpha=0.5, linewidth=3)

axs[1].set_xlabel(r"$\sigma_{\rm HI}$ (km/s)")
axs[0].set_ylabel(r"$\sigma_{\rm CO}$ (km/s)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_cloudtype_nofit.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_vs_H2_w_cloudtype_nofit.pdf"))
plt.close()


# Make histograms of each distribution and the cloud types
from astropy.stats import histogram

onecolumn_figure()

vals, bins = histogram(tab['sigma_HI'][good_pts] / 1000.,
                       bins='knuth')
vals_a = histogram(tab['sigma_HI'][good_pts & (tab['cloud_type'] == "A")] / 1000.,
                   bins=bins)[0]
vals_b = histogram(tab['sigma_HI'][good_pts & (tab['cloud_type'] == "B")] / 1000.,
                   bins=bins)[0]
vals_c = histogram(tab['sigma_HI'][good_pts & (tab['cloud_type'] == "C")] / 1000.,
                   bins=bins)[0]

plt.plot((bins[1:] + bins[:-1]) * 0.5, vals, 'k', linewidth=3,
         drawstyle='steps-mid')

plt.grid()
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_HI_hist.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_hist.pdf"))
plt.close()

plt.plot((bins[1:] + bins[:-1]) * 0.5, vals, 'k', linewidth=3,
         drawstyle='steps-mid')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_a, linewidth=1.5,
         drawstyle='steps-mid', label='A')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_b, linewidth=1.5,
         drawstyle='steps-mid', label='B')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_c, linewidth=1.5,
         drawstyle='steps-mid', label='C')

plt.grid()
plt.legend(frameon=True)
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_HI_hist_wcloudtype.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_hist_wcloudtype.pdf"))
plt.close()

# Make a second version in density

vals, bins = histogram(tab['sigma_HI'][good_pts] / 1000.,
                       bins='knuth', density=True)
vals_a = histogram(tab['sigma_HI'][good_pts & (tab['cloud_type'] == "A")] / 1000.,
                   bins=bins, density=True)[0]
vals_b = histogram(tab['sigma_HI'][good_pts & (tab['cloud_type'] == "B")] / 1000.,
                   bins=bins, density=True)[0]
vals_c = histogram(tab['sigma_HI'][good_pts & (tab['cloud_type'] == "C")] / 1000.,
                   bins=bins, density=True)[0]

plt.plot((bins[1:] + bins[:-1]) * 0.5, vals, 'k', linewidth=3,
         drawstyle='steps-mid')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_a, linewidth=1.5,
         drawstyle='steps-mid', label='A')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_b, linewidth=1.5,
         drawstyle='steps-mid', label='B')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_c, linewidth=1.5,
         drawstyle='steps-mid', label='C')

plt.grid()
plt.legend(frameon=True)
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_HI_hist_density_wcloudtype.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_hist_density_wcloudtype.pdf"))
plt.close()
# Now make a CDF plot

plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals) / np.sum(vals),
         'k', linewidth=3,
         drawstyle='steps-mid')
plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals_a) / np.sum(vals_a),
         linewidth=1.5,
         drawstyle='steps-mid', label='A')
plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals_b) / np.sum(vals_b),
         linewidth=1.5,
         drawstyle='steps-mid', label='B')
plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals_c) / np.sum(vals_c),
         linewidth=1.5,
         drawstyle='steps-mid', label='C')

plt.grid()
plt.legend(frameon=True)
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_HI_cdf_wcloudtype.png"))
plt.savefig(osjoin(fig_path, "sigma_HI_cdf_wcloudtype.pdf"))
plt.close()

# And CO

vals, bins = histogram(tab['sigma_CO'][good_pts] / 1000.,
                       bins='knuth')
vals_a = histogram(tab['sigma_CO'][good_pts & (tab['cloud_type'] == "A")] / 1000.,
                   bins=bins)[0]
vals_b = histogram(tab['sigma_CO'][good_pts & (tab['cloud_type'] == "B")] / 1000.,
                   bins=bins)[0]
vals_c = histogram(tab['sigma_CO'][good_pts & (tab['cloud_type'] == "C")] / 1000.,
                   bins=bins)[0]

plt.plot((bins[1:] + bins[:-1]) * 0.5, vals, 'k', linewidth=3,
         drawstyle='steps-mid')
plt.grid()
plt.xlabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_CO_hist.png"))
plt.savefig(osjoin(fig_path, "sigma_CO_hist.pdf"))
plt.close()

plt.plot((bins[1:] + bins[:-1]) * 0.5, vals, 'k', linewidth=3,
         drawstyle='steps-mid')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_a, linewidth=1.5,
         drawstyle='steps-mid', label='A')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_b, linewidth=1.5,
         drawstyle='steps-mid', label='B')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_c, linewidth=1.5,
         drawstyle='steps-mid', label='C')

plt.grid()
plt.legend(frameon=True)
plt.xlabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_CO_hist_wcloudtype.png"))
plt.savefig(osjoin(fig_path, "sigma_CO_hist_wcloudtype.pdf"))
plt.close()

# And density

vals, bins = histogram(tab['sigma_CO'][good_pts] / 1000.,
                       bins='knuth', density=True)
vals_a = histogram(tab['sigma_CO'][good_pts & (tab['cloud_type'] == "A")] / 1000.,
                   bins=bins, density=True)[0]
vals_b = histogram(tab['sigma_CO'][good_pts & (tab['cloud_type'] == "B")] / 1000.,
                   bins=bins, density=True)[0]
vals_c = histogram(tab['sigma_CO'][good_pts & (tab['cloud_type'] == "C")] / 1000.,
                   bins=bins, density=True)[0]

plt.plot((bins[1:] + bins[:-1]) * 0.5, vals, 'k', linewidth=3,
         drawstyle='steps-mid')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_a, linewidth=1.5,
         drawstyle='steps-mid', label='A')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_b, linewidth=1.5,
         drawstyle='steps-mid', label='B')
plt.plot((bins[1:] + bins[:-1]) * 0.5, vals_c, linewidth=1.5,
         drawstyle='steps-mid', label='C')

plt.grid()
plt.legend(frameon=True)
plt.xlabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_CO_hist_density_wcloudtype.png"))
plt.savefig(osjoin(fig_path, "sigma_CO_hist_density_wcloudtype.pdf"))
plt.close()

# CO CDF

plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals) / np.sum(vals),
         'k', linewidth=3,
         drawstyle='steps-mid')
plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals_a) / np.sum(vals_a),
         linewidth=1.5,
         drawstyle='steps-mid', label='A')
plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals_b) / np.sum(vals_b),
         linewidth=1.5,
         drawstyle='steps-mid', label='B')
plt.plot((bins[1:] + bins[:-1]) * 0.5,
         np.cumsum(vals_c) / np.sum(vals_c),
         linewidth=1.5,
         drawstyle='steps-mid', label='C')

plt.grid()
plt.legend(frameon=True)
plt.xlabel(r"$\sigma_{\rm CO}$ (km/s)")
plt.tight_layout()
plt.savefig(osjoin(fig_path, "sigma_CO_cdf_wcloudtype.png"))
plt.savefig(osjoin(fig_path, "sigma_CO_cdf_wcloudtype.pdf"))
plt.close()
