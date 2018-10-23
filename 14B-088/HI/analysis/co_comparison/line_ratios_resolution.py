
'''
Explore the line width ratios with different spatial resolutions.
'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import seaborn as sb
import emcee
from corner import hist2d

from plotting_styles import default_figure, onecolumn_figure, onecolumn_Npanel_figure, twocolumn_figure

from paths import (iram_co21_14B088_data_path, fourteenB_HI_data_wGBT_path, allfigs_path)


# Load the radial-stacked tables in
hi_tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_radial_500pc.csv"))
hi_tab_2beam = Table.read(fourteenB_HI_data_wGBT_path("smooth_2beam/tables/hi_hwhm_totalprof_fits_38arcsec_radial_500pc.csv"))
hi_tab_5beam = Table.read(fourteenB_HI_data_wGBT_path("smooth_5beam/tables/hi_hwhm_totalprof_fits_95arcsec_radial_500pc.csv"))

co_tab = Table.read(iram_co21_14B088_data_path("tables/co_hwhm_totalprof_fits_radial_500pc.csv"))
co_tab_2beam = Table.read(iram_co21_14B088_data_path("smooth_2beam/tables/co_hwhm_totalprof_fits_38arcsec_radial_500pc.csv"))
co_tab_5beam = Table.read(iram_co21_14B088_data_path("smooth_5beam/tables/co_hwhm_totalprof_fits_95arcsec_radial_500pc.csv"))

bin_centers = np.arange(0., 6.6, 0.5) + 0.25

onecolumn_Npanel_figure(N=2)

fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)

ax0 = ax[0]
ax1 = ax[1]

ax0.errorbar(bin_centers, hi_tab['peaksub_sigma'],
             yerr=hi_tab['peaksub_sigma_low_lim'], label='80 pc (19")',
             drawstyle='steps-mid')
ax0.errorbar(bin_centers, hi_tab_2beam['peaksub_sigma'],
             yerr=hi_tab_2beam['peaksub_sigma_low_lim'], label='160 pc (38")',
             drawstyle='steps-mid', linestyle='--')
ax0.errorbar(bin_centers, hi_tab_5beam['peaksub_sigma'],
             yerr=hi_tab_5beam['peaksub_sigma_low_lim'], label='380 pc (95")',
             drawstyle='steps-mid', linestyle=':')
ax0.set_ylabel(r"$\sigma_{\rm HI}$ (km/s)")
ax0.grid()
ax0.legend(frameon=True)
ax0.set_ylim([0, 16])

ax1.errorbar(bin_centers, co_tab['peaksub_sigma'],
             yerr=co_tab['peaksub_sigma_low_lim'], label='19"',
             drawstyle='steps-mid')
ax1.errorbar(bin_centers, co_tab_2beam['peaksub_sigma'],
             yerr=co_tab_2beam['peaksub_sigma_low_lim'], label='38"',
             drawstyle='steps-mid', linestyle='--')
ax1.errorbar(bin_centers, co_tab_5beam['peaksub_sigma'],
             yerr=co_tab_5beam['peaksub_sigma_low_lim'], label='95"',
             drawstyle='steps-mid', linestyle=':')
ax1.set_ylabel(r"$\sigma_{\rm CO}$ (km/s)")

ax1.grid()

ax1.set_xlabel("Radius (kpc)")

plt.tight_layout()

plt.savefig(allfigs_path("stacked_profiles/total_profile_radial_widths_HI_CO21_varyres.pdf"))
plt.savefig(allfigs_path("stacked_profiles/total_profile_radial_widths_HI_CO21_varyres.png"))
plt.close()

ratio = co_tab['peaksub_sigma'] / hi_tab['peaksub_sigma']
ratio_err = ratio * np.sqrt((co_tab['peaksub_sigma_low_lim'] / co_tab['peaksub_sigma'])**2 +
                            (hi_tab['peaksub_sigma_low_lim'] / hi_tab['peaksub_sigma'])**2)

ratio_2beam = co_tab_2beam['peaksub_sigma'] / hi_tab_2beam['peaksub_sigma']
ratio_2beam_err = ratio_2beam * \
    np.sqrt((co_tab_2beam['peaksub_sigma_low_lim'] / co_tab_2beam['peaksub_sigma'])**2 +
            (hi_tab_2beam['peaksub_sigma_low_lim'] / hi_tab_2beam['peaksub_sigma'])**2)

ratio_5beam = co_tab_5beam['peaksub_sigma'] / hi_tab_5beam['peaksub_sigma']
ratio_5beam_err = ratio_5beam * \
    np.sqrt((co_tab_5beam['peaksub_sigma_low_lim'] / co_tab_5beam['peaksub_sigma'])**2 +
            (hi_tab_5beam['peaksub_sigma_low_lim'] / hi_tab_5beam['peaksub_sigma'])**2)

onecolumn_figure()

plt.errorbar(bin_centers, ratio,
             yerr=ratio_err, label='19"',
             drawstyle='steps-mid')
plt.errorbar(bin_centers, ratio_2beam,
             yerr=ratio_2beam_err, label='38"',
             drawstyle='steps-mid')
plt.errorbar(bin_centers, ratio_5beam,
             yerr=ratio_5beam_err, label='95"',
             drawstyle='steps-mid')
plt.ylabel(r"$\sigma_{\rm CO}$ / $\sigma_{\rm HI}$")
plt.xlabel("Radius (kpc)")
plt.grid()

plt.tight_layout()

plt.savefig(allfigs_path("stacked_profiles/total_profile_radial_widthratios_HI_CO21_varyres.pdf"))
plt.savefig(allfigs_path("stacked_profiles/total_profile_radial_widthratios_HI_CO21_varyres.png"))
plt.close()

# Load the LOS tables

hi_co_tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_co_gaussfit_column_densities_perpix.fits"))
hi_co_tab_2beam = Table.read(fourteenB_HI_data_wGBT_path("smooth_2beam/tables/hi_co_gaussfit_column_densities_perpix_38arcsec.fits"))
hi_co_tab_5beam = Table.read(fourteenB_HI_data_wGBT_path("smooth_5beam/tables/hi_co_gaussfit_column_densities_perpix_95arcsec.fits"))

good_pts = np.logical_and(~hi_co_tab['multicomp_flag_HI'],
                          ~hi_co_tab['multicomp_flag_CO'])
good_pts = np.logical_and(good_pts,
                          hi_co_tab["sigma_HI"] > 3800)

good_pts_2beam = np.logical_and(~hi_co_tab_2beam['multicomp_flag_HI'],
                                ~hi_co_tab_2beam['multicomp_flag_CO'])
# good_pts_2beam = np.logical_and(good_pts_2beam,
#                                 hi_co_tab_2beam["sigma_HI"] > 3800)

good_pts_5beam = np.logical_and(~hi_co_tab_5beam['multicomp_flag_HI'],
                                ~hi_co_tab_5beam['multicomp_flag_CO'])
# good_pts_5beam = np.logical_and(good_pts_5beam,
#                                 hi_co_tab_5beam["sigma_HI"] > 3800)

# What are the median line widths
print("19'' 15th, 50th, 85th line widths. HI: {0}; CO: {1}"
      .format(np.percentile(hi_co_tab['sigma_HI'][good_pts], [15, 50, 85]),
              np.percentile(hi_co_tab['sigma_CO'][good_pts], [15, 50, 85])))
print("38'' 15th, 50th, 85th line widths. HI: {0}; CO: {1}"
      .format(np.percentile(hi_co_tab_2beam['sigma_HI'][good_pts_2beam], [15, 50, 85]),
              np.percentile(hi_co_tab_2beam['sigma_CO'][good_pts_2beam], [15, 50, 85])))
print("95'' 15th, 50th, 85th line widths. HI: {0}; CO: {1}"
      .format(np.percentile(hi_co_tab_5beam['sigma_HI'][good_pts_5beam], [15, 50, 85]),
              np.percentile(hi_co_tab_5beam['sigma_CO'][good_pts_5beam], [15, 50, 85])))
# 19'' 15th, 50th, 85th line widths.
# HI: [ 6151.29135369  7429.32916764  9175.82252398];
# CO: [ 3289.7303243   4288.90275791  5801.99530487]
# 38'' 15th, 50th, 85th line widths.
# HI: [  7179.47168223   8358.83189562  10151.71567788];
# CO: [ 3935.75509445  4999.51314003  6359.43123789]
# 95'' 15th, 50th, 85th line widths.
# HI: [  9071.86256852  11022.19067997  13678.21229086];
# CO: [ 5795.67825932  7419.73937683  9789.86425635]


twocolumn_figure()

fig, axes = plt.subplots(2, 3, sharex=True, sharey='row')

sb.violinplot(x=(hi_co_tab['sigma_HI'][good_pts] / 1000.).tolist(),
              ax=axes[0, 0], orient='v')
sb.violinplot(x=(hi_co_tab_2beam['sigma_HI'][good_pts_2beam] / 1000.).tolist(),
              ax=axes[0, 1], orient='v')
sb.violinplot(x=(hi_co_tab_5beam['sigma_HI'][good_pts_5beam] / 1000.).tolist(),
              ax=axes[0, 2], orient='v')

sb.violinplot(x=(hi_co_tab['sigma_CO'][good_pts] / 1000.).tolist(),
              ax=axes[1, 0], orient='v')
sb.violinplot(x=(hi_co_tab_2beam['sigma_CO'][good_pts_2beam] / 1000.).tolist(),
              ax=axes[1, 1], orient='v')
sb.violinplot(x=(hi_co_tab_5beam['sigma_CO'][good_pts_5beam] / 1000.).tolist(),
              ax=axes[1, 2], orient='v')

axes[0, 0].set_ylabel(r"$\sigma_{\rm HI}$ (km/s)")
axes[1, 0].set_ylabel(r"$\sigma_{\rm CO}$ (km/s)")

axes[1, 0].set_xlabel(r"$19''$")
axes[1, 1].set_xlabel(r"$38''$")
axes[1, 2].set_xlabel(r"$95''$")

for ax in axes.ravel():
    ax.grid()

fig.tight_layout()

fig.savefig(allfigs_path("co_vs_hi/sigma_HI_CO_fit_varyres.pdf"))
fig.savefig(allfigs_path("co_vs_hi/sigma_HI_CO_fit_varyres.png"))
plt.close()

print(np.median(hi_co_tab['sigma_HI'][good_pts]),
      np.median(hi_co_tab['sigma_CO'][good_pts]),
      np.median(hi_co_tab['sigma_CO'][good_pts]) /
      np.median(hi_co_tab['sigma_HI'][good_pts]))

print(np.median(hi_co_tab_2beam['sigma_HI'][good_pts_2beam]),
      np.median(hi_co_tab_2beam['sigma_CO'][good_pts_2beam]),
      np.median(hi_co_tab_2beam['sigma_CO'][good_pts_2beam]) /
      np.median(hi_co_tab_2beam['sigma_HI'][good_pts_2beam]))

print(np.median(hi_co_tab_5beam['sigma_HI'][good_pts_5beam]),
      np.median(hi_co_tab_5beam['sigma_CO'][good_pts_5beam]),
      np.median(hi_co_tab_5beam['sigma_CO'][good_pts_5beam]) /
      np.median(hi_co_tab_5beam['sigma_HI'][good_pts_5beam]))

# Find the overlapping set of points between all.
overlap_mask = np.zeros(len(hi_co_tab["mean_HI"]), dtype=bool)

for y, x in zip(hi_co_tab_5beam['ypts'], hi_co_tab_5beam['xpts']):
    yidx = np.where(hi_co_tab['ypts'] == y)[0]
    xidx = np.where(hi_co_tab['xpts'] == x)[0]

    # Look for overlaps
    match = list(set(yidx) & set(xidx))

    if len(match) == 0:
        continue
    if len(match) > 1:
        raise ValueError("Multiple matches? Not possible!")

    overlap_mask[match[0]] = True

overlap_mask_5 = np.zeros(len(hi_co_tab_2beam["mean_HI"]), dtype=bool)

for y, x in zip(hi_co_tab_5beam['ypts'], hi_co_tab_5beam['xpts']):
    yidx = np.where(hi_co_tab_2beam['ypts'] == y)[0]
    xidx = np.where(hi_co_tab_2beam['xpts'] == x)[0]

    # Look for overlaps
    match = list(set(yidx) & set(xidx))

    if len(match) == 0:
        continue
    if len(match) > 1:
        raise ValueError("Multiple matches? Not possible!")

    overlap_mask_5[match[0]] = True

overlap_mask_2 = np.zeros(len(hi_co_tab["mean_HI"]), dtype=bool)

for y, x in zip(hi_co_tab_2beam['ypts'], hi_co_tab_2beam['xpts']):
    yidx = np.where(hi_co_tab['ypts'] == y)[0]
    xidx = np.where(hi_co_tab['xpts'] == x)[0]

    # Look for overlaps
    match = list(set(yidx) & set(xidx))

    if len(match) == 0:
        continue
    if len(match) > 1:
        raise ValueError("Multiple matches? Not possible!")

    overlap_mask_2[match[0]] = True

good_pts_orig = np.logical_and(good_pts, overlap_mask)

good_pts_2 = np.logical_and(good_pts[overlap_mask_2], overlap_mask_5)

good_pts_5 = good_pts[overlap_mask]

good_pts_5_only = np.logical_and(~hi_co_tab_5beam['multicomp_flag_HI'],
                                 ~hi_co_tab_5beam['multicomp_flag_CO'])

good_pts_2_only = np.logical_and(~hi_co_tab_2beam['multicomp_flag_HI'],
                                 ~hi_co_tab_2beam['multicomp_flag_CO'])


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


# Fit the smoothed relations

# params, cis, sampler = \
#     bayes_linear(hi_co_tab_2beam['sigma_HI'][good_pts_2beam],
#                  hi_co_tab_2beam['sigma_CO'][good_pts_2beam],
#                  hi_co_tab_2beam['sigma_stderr_HI'][good_pts_2beam],
#                  hi_co_tab_2beam['sigma_stderr_CO'][good_pts_2beam],
#                  nBurn=500, nSample=2000, nThin=2)

# Just fit ratio (intercept is 0)
pars_ratio, pars_ratio_ci, sampler_ratio = \
    bayes_linear(hi_co_tab_2beam['sigma_HI'][good_pts_2beam],
                 hi_co_tab_2beam['sigma_CO'][good_pts_2beam],
                 hi_co_tab_2beam['sigma_stderr_HI'][good_pts_2beam],
                 hi_co_tab_2beam['sigma_stderr_CO'][good_pts_2beam],
                 nBurn=500, nSample=5000, nThin=1,
                 fix_intercept=True)

slope_ratio = pars_ratio[0]
slope_ratio_ci = pars_ratio_ci[0]

add_stddev_ratio = pars_ratio[1]
add_stddev_ratio_ci = pars_ratio_ci[1]

onecolumn_figure()
hist2d(hi_co_tab_2beam['sigma_HI'][good_pts_2beam] / 1000.,
       hi_co_tab_2beam['sigma_CO'][good_pts_2beam] / 1000., bins=10,
       data_kwargs={"alpha": 0.5})
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")

# slope = params[0]
# inter = params[1] / 1000.
# slope_ci = cis[0]
# inter_cis = cis[1] / 1000.
# plt.plot([4, 12], [4. * slope + inter, 12. * slope + inter])
# plt.fill_between([4, 12], [4. * slope_ci[0] + inter_cis[0],
#                            12. * slope_ci[0] + inter_cis[0]],
#                  [4. * slope_ci[1] + inter_cis[1],
#                   12. * slope_ci[1] + inter_cis[1]],
#                  facecolor=sb.color_palette()[0],
#                  alpha=0.5)
plt.plot([4, 12], [4. * slope_ratio, 12. * slope_ratio],
         '--', color=sb.color_palette()[1], linewidth=3)

plt.plot([4, 12], [4, 12], '-.', linewidth=3, alpha=0.8,
         color=sb.color_palette()[2])

plt.tight_layout()

plt.savefig(allfigs_path("co_vs_hi/sigma_HI_vs_H2_w_fit_38arcsec.png"))
plt.savefig(allfigs_path("co_vs_hi/sigma_HI_vs_H2_w_fit_38arcsec.pdf"))
plt.close()

# print("2-beam Slope: {0} {1}".format(slope, slope_ci))
# print("2-beam Intercept: {0} {1}".format(inter, inter_cis))

print("2-beam Ratio Slope: {0} {1}".format(slope_ratio, slope_ratio_ci))
# 2-beam Ratio Slope: 0.566189501877 [ 0.56529157  0.56708893]

print("2-beam Added stddev: {0} {1}".format(add_stddev_ratio, add_stddev_ratio_ci))
# Added stddev: 302.786781711 [ 294.64176568  310.5069717 ]


# params_5, cis_5, sampler_5 = \
#     bayes_linear(hi_co_tab_5beam['sigma_HI'][good_pts_5beam],
#                  hi_co_tab_5beam['sigma_CO'][good_pts_5beam],
#                  hi_co_tab_5beam['sigma_stderr_HI'][good_pts_5beam],
#                  hi_co_tab_5beam['sigma_stderr_CO'][good_pts_5beam],
#                  nBurn=500, nSample=2000, nThin=2)

# Just fit ratio (intercept is 0)
pars_ratio_5, pars_ratio_ci_5, sampler_ratio_5 = \
    bayes_linear(hi_co_tab_5beam['sigma_HI'][good_pts_5beam],
                 hi_co_tab_5beam['sigma_CO'][good_pts_5beam],
                 hi_co_tab_5beam['sigma_stderr_HI'][good_pts_5beam],
                 hi_co_tab_5beam['sigma_stderr_CO'][good_pts_5beam],
                 nBurn=500, nSample=5000, nThin=1,
                 fix_intercept=True)

slope_ratio_5 = pars_ratio_5[0]
slope_ratio_ci_5 = pars_ratio_ci_5[0]

add_stddev_ratio_5 = pars_ratio_5[1]
add_stddev_ratio_ci_5 = pars_ratio_ci_5[1]

onecolumn_figure()
hist2d(hi_co_tab_5beam['sigma_HI'][good_pts_5beam] / 1000.,
       hi_co_tab_5beam['sigma_CO'][good_pts_5beam] / 1000., bins=10,
       data_kwargs={"alpha": 0.5})
plt.xlabel(r"$\sigma_{\rm HI}$ (km/s)")
plt.ylabel(r"$\sigma_{\rm CO}$ (km/s)")

# slope_5 = params_5[0]
# inter_5 = params_5[1] / 1000.
# slope_ci_5 = cis_5[0]
# inter_cis_5 = cis_5[1] / 1000.
# plt.plot([4, 17], [4. * slope_5 + inter_5, 17. * slope_5 + inter_5])
# plt.fill_between([4, 17], [4. * slope_ci_5[0] + inter_cis_5[0],
#                            17. * slope_ci_5[0] + inter_cis_5[0]],
#                  [4. * slope_ci_5[1] + inter_cis_5[1],
#                   17. * slope_ci_5[1] + inter_cis_5[1]],
#                  facecolor=sb.color_palette()[0],
#                  alpha=0.5)
plt.plot([4, 17], [4. * slope_ratio_5, 17. * slope_ratio_5],
         '--', color=sb.color_palette()[1], linewidth=3)

plt.plot([4, 17], [4, 17], '-.', linewidth=3, alpha=0.8,
         color=sb.color_palette()[2])

plt.tight_layout()
plt.savefig(allfigs_path("co_vs_hi/sigma_HI_vs_H2_w_fit_95arcsec.png"))
plt.savefig(allfigs_path("co_vs_hi/sigma_HI_vs_H2_w_fit_95arcsec.pdf"))
plt.close()

# print("5-beam Slope: {0} {1}".format(slope_5, slope_ci_5))
# print("5-beam Intercept: {0} {1}".format(inter_5, inter_cis_5))

print("5-beam Ratio Slope: {0} {1}".format(slope_ratio_5, slope_ratio_ci_5))
# 5-beam Ratio Slope: 0.632291388013 [ 0.63100664  0.63355342]

print("5-beam Added stddev: {0} {1}".format(add_stddev_ratio_5,
                                            add_stddev_ratio_ci_5))
# Added stddev: 23.4024470066 [ 11.42358035  38.02678363]
