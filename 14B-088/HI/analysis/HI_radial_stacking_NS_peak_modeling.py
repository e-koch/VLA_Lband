
'''
Compare properties of the peak velocity stacks between the N/S halves
relative to the full stack.
'''
import astropy.units as u
from spectral_cube import SpectralCube, OneDSpectrum
import numpy as np
import matplotlib.pyplot as p
from pandas import read_csv, DataFrame
from astropy.io import fits
from os.path import join as osjoin
import os
import seaborn as sb

from cube_analysis.spectral_stacking_models import (fit_hwhm,
                                                    find_linewing_asymm)

from paths import (allfigs_path, alltables_path, fourteenB_HI_data_path,
                   fourteenB_HI_data_wGBT_path)
from plotting_styles import *

'''
Create profiles of HI after subtracting velocity surfaces.
'''

noSD_stackpath = lambda x: osjoin(fourteenB_HI_data_path("", no_check=True), "stacked_spectra", x)
comb_stackpath = lambda x: osjoin(fourteenB_HI_data_wGBT_path("", no_check=True), "stacked_spectra", x)


figure_folder = "stacked_profiles"
if not os.path.exists(allfigs_path(figure_folder)):
    os.mkdir(allfigs_path(figure_folder))

# Here are the number of pixels per bin for the different masks
npix_beam = 41.

num_pix_total = 818998
num_pix_feath_total = 965943

num_pix = \
    np.array([118., 343., 574., 806., 1053., 1270., 1500.,
             1729., 1965., 2185., 2420., 2676., 2895., 3117.,
             3348., 3563., 3813., 4045., 4285., 4510., 4753.,
             4949., 5185., 5424., 5678., 5898., 6126., 6345.,
             6586., 6820., 7018., 7284., 7532., 7754., 7963.,
             8188., 8434., 8656., 8892., 9125., 9379., 9607.,
             9794., 10053., 10268., 10507., 10742., 10992., 11218.,
             11435., 11646., 11887., 12121., 12367., 12620., 12835.,
             13045., 13272., 13499., 13735., 14001., 14230., 14425.,
             14686., 14897., 15119., 15343., 15599., 15852., 16060.,
             16302., 16480., 16771., 16948., 17206., 17471., 17685.,
             17917., 18114., 18358.])

n_numpix = \
    np.array([59., 171., 288., 402., 528., 637., 745., 866.,
              985., 1091., 1204., 1340., 1450., 1560., 1679., 1771.,
              1911., 2015., 2149., 2257., 2371., 2478., 2593., 2716.,
              2828., 2954., 3065., 3176., 3297., 3405., 3509., 3628.,
              3779., 3875., 3986., 4091., 4221., 4336., 4436., 4548.,
              4690., 4823., 4894., 5027., 5130., 5252., 5366., 5499.,
              5621., 5702., 5825., 5942., 6072., 6169., 6310., 6421.,
              6526., 6640., 6745., 6866., 6991., 7125., 7225., 7336.,
              7449., 7553., 7687., 7775., 7940., 8026., 8155., 8246.,
              8387., 8469., 8582., 8746., 8852., 8972., 9047., 9174.])

s_numpix = \
    np.array([59., 172., 286., 404., 525., 633., 755., 863.,
              980., 1094., 1216., 1336., 1445., 1557., 1669., 1792.,
              1902., 2030., 2136., 2253., 2382., 2471., 2592., 2708.,
              2850., 2944., 3061., 3169., 3289., 3415., 3509., 3656.,
              3753., 3879., 3977., 4097., 4213., 4320., 4456., 4577.,
              4689., 4784., 4900., 5026., 5138., 5255., 5376., 5493.,
              5597., 5733., 5821., 5945., 6049., 6198., 6310., 6414.,
              6519., 6632., 6754., 6869., 7010., 7105., 7200., 7350.,
              7448., 7566., 7656., 7824., 7912., 8034., 8147., 8234.,
              8384., 8479., 8624., 8725., 8833., 8945., 9067., 9184.])

sigma_noise = 2.8  # K

dr = 100 * u.pc
max_radius = (8.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)

nbins = np.int(np.floor(max_radius / dr))
inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)
bin_cents = (outeredge - dr / 2.).to(u.kpc).value

total_peakvel = OneDSpectrum.from_hdu(fits.open(noSD_stackpath("peakvel_stacked.fits")))
total_peakvel_wGBT = OneDSpectrum.from_hdu(fits.open(comb_stackpath("peakvel_stacked.fits")))

peakvel_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring)))
peakvel_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring)))
peakvel_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring)))

peakvel_stack_feath_n = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring)))
peakvel_stack_feath_s = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring)))
peakvel_stack_feath = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring)))

total_peakvel_n = peakvel_stack_n.sum(axis=(1, 2))
total_peakvel_s = peakvel_stack_s.sum(axis=(1, 2))

total_peakvel_feath_n = peakvel_stack_feath_n.sum(axis=(1, 2))
total_peakvel_feath_s = peakvel_stack_feath_s.sum(axis=(1, 2))

# Saved fits from HI_stacking_modeling
hi_radial_fits = read_csv(fourteenB_HI_data_path("tables/hi_gaussian_hwhm_totalprof_fits_radial.csv"),
                          index_col=0)

# First plots comparing the total to N and S
spectra = [total_peakvel, total_peakvel_wGBT,
           total_peakvel_n, total_peakvel_s,
           total_peakvel_feath_n, total_peakvel_feath_s]
labels = ["Peak Sub.", "Feath. Peak Sub.",
          "Peak Sub. N", 'Peak Sub. S',
          "Feath. Peak Sub. N", "Feath. Peak Sub. S"]
file_labels = ["peaksub", "peaksub_feath",
               "peaksub_n", "peaksub_s",
               "peaksub_feath_n", "peaksub_feath_s"]
hi_fit_hwhm_vals = {}
hwhm_models = dict.fromkeys(file_labels)

onecolumn_figure()

for spectrum, label, file_label in zip(spectra, labels, file_labels):

    norm_intens = (spectrum / np.nanmax(spectrum)).value

    vels = spectrum.spectral_axis.value / 1000.

    # HWHM fitting
    parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_HI_hwhm = \
        fit_hwhm(vels, norm_intens,
                 sigma_noise=sigma_noise / np.nanmax(spectrum.value),
                 nbeams=(num_pix_feath_total / npix_beam if "feath" in
                         file_label else num_pix_total / npix_beam),
                 niters=100, interp_factor=2.)

    hi_fit_hwhm_vals[label + " Params"] = parvals_hwhm
    hi_fit_hwhm_vals[label + " Lower Limit"] = np.abs(parerrs_hwhm[0])
    hi_fit_hwhm_vals[label + " Upper Limit"] = np.abs(parerrs_hwhm[1])
    hwhm_models[file_label] = g_HI_hwhm

    # Note that the statistical errors on the mean are too small.

    # Skip the non-N or S stacks. Already saved in HI_stacking_modeling.py
    if "_n" not in file_label and "_s" not in file_label:
        continue

    # Now plot the HWHM model.
    ax1 = p.subplot2grid((5, 1), (0, 0), rowspan=4)
    p.plot(vels, norm_intens, '-', drawstyle='steps-mid')
    p.plot(vels, g_HI_hwhm(vels), ':', label="Fit")
    p.fill_between(vels, g_HI_hwhm(vels), norm_intens, facecolor='gray',
                   alpha=0.5)
    p.ylabel("HI Normalized Intensity")
    p.xticks([])
    p.xlim([-50, 50])
    p.legend(frameon=True)
    p.ylim([-0.1, 1.1])
    p.grid()

    ax2 = p.subplot2grid((5, 1), (4, 0))
    p.plot(vels, norm_intens - g_HI_hwhm(vels), '-', drawstyle='steps-mid')
    p.hlines(0.0, -100, 100, color='k')
    p.grid()
    # p.ylim([-0.015, 0.015])
    p.xlim([-50, 50])
    p.xlabel("Velocity (km/s)")
    p.ylabel("Residuals")

    p.tight_layout()
    p.subplots_adjust(hspace=0)

    filename = "total_profile_corrected_velocity_{}_hi_hwhm_fit".format(file_label)
    p.savefig(allfigs_path(osjoin(figure_folder, filename + ".pdf")))
    p.savefig(allfigs_path(osjoin(figure_folder, filename + ".png")))
    p.close()
    # p.draw()
    # raw_input("Next plot?")
    # p.clf()

spectra = [total_peakvel, total_peakvel_wGBT,
           total_peakvel_n, total_peakvel_s,
           total_peakvel_feath_n, total_peakvel_feath_s]

# Plot comparison between total profiles
onecolumn_figure()

pos_vel_mask = total_peakvel_wGBT.spectral_axis.to(u.km / u.s).value > 0.
neg_vel_mask = total_peakvel_wGBT.spectral_axis.to(u.km / u.s).value < 0.

p.plot(total_peakvel_feath_n.spectral_axis.to(u.km / u.s).value,
       (total_peakvel_feath_n / total_peakvel_feath_n.max()).value,
       '--', drawstyle='steps-mid', label='North',
           color=sb.color_palette()[2])
p.plot(total_peakvel_feath_s.spectral_axis.to(u.km / u.s).value,
       (total_peakvel_feath_s / total_peakvel_feath_s.max()).value,
       '-.', drawstyle='steps-mid', label='South',
       color=sb.color_palette()[0])

# Fill the differences in the line wings
p.fill_between(total_peakvel_feath_s.spectral_axis.to(u.km / u.s).value[pos_vel_mask],
               (total_peakvel_feath_s / total_peakvel_feath_s.max()).value[pos_vel_mask],
               (total_peakvel_feath_n / total_peakvel_feath_n.max()).value[pos_vel_mask],
               color=sb.color_palette()[2], alpha=0.8)
p.fill_between(total_peakvel_feath_s.spectral_axis.to(u.km / u.s).value[neg_vel_mask],
               (total_peakvel_feath_s / total_peakvel_feath_s.max()).value[neg_vel_mask],
               (total_peakvel_feath_n / total_peakvel_feath_n.max()).value[neg_vel_mask],
               color=sb.color_palette()[0], alpha=0.8)
p.plot(total_peakvel_wGBT.spectral_axis.to(u.km / u.s).value,
       (total_peakvel_wGBT / total_peakvel_wGBT.max()).value,
       '-', drawstyle='steps-mid', label='Total',
       color='k', alpha=0.7)

# p.text(-40, 0.88, "VLA+GBT",
#          bbox={"boxstyle": "square", "facecolor": "w"})
p.ylim([-0.02, 1.03])
p.xlim([-50, 50])
p.grid()
p.legend(frameon=True, loc='upper right')
p.xlabel("Velocity (km/s)")
p.ylabel("Normalized Total Intensity")

p.tight_layout()

p.savefig(allfigs_path(osjoin(figure_folder, "total_profile_NS_comparison.pdf")))
p.savefig(allfigs_path(osjoin(figure_folder, "total_profile_NS_comparison.png")))
p.close()

onecolumn_Npanel_figure(N=2.3)
fig, ax = p.subplots(3, 1, sharey=False, sharex=True)

ax[2].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_f_wings'],
               yerr=[hi_radial_fits["peaksub_feath_f_wings_low_lim"],
                     hi_radial_fits["peaksub_feath_f_wings_up_lim"]],
               linestyle='-', drawstyle='steps-mid',
               color='k', alpha=0.7)
ax[2].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_n_f_wings'],
               yerr=[hi_radial_fits["peaksub_feath_n_f_wings_low_lim"],
                     hi_radial_fits["peaksub_feath_n_f_wings_up_lim"]],
               linestyle='--', drawstyle='steps-mid',
               color=sb.color_palette()[2])
ax[2].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_s_f_wings'],
               yerr=[hi_radial_fits["peaksub_feath_s_f_wings_low_lim"],
                     hi_radial_fits["peaksub_feath_s_f_wings_up_lim"]],
               linestyle='-.', drawstyle='steps-mid',
               color=sb.color_palette()[0])

ax[2].set_xlim([0, 8])
ax[2].set_ylim([0.13, 0.38])
ax[2].grid()
ax[2].set_ylabel(r"$f_{\rm wings}$")
ax[2].set_xlabel("Radius (kpc)")
ax[2].fill_between([0, 0.5], 0.13, 0.38, facecolor='gray', alpha=0.5,
                   zorder=-1)

# kappa
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_kappa'],
               yerr=[hi_radial_fits["peaksub_feath_kappa_low_lim"],
                     hi_radial_fits["peaksub_feath_kappa_up_lim"]],
               linestyle='-', drawstyle='steps-mid',
               color='k', alpha=0.7)
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_n_kappa'],
               yerr=[hi_radial_fits["peaksub_feath_n_kappa_low_lim"],
                     hi_radial_fits["peaksub_feath_n_kappa_up_lim"]],
               linestyle='--', drawstyle='steps-mid',
               color=sb.color_palette()[2])
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_s_kappa'],
               yerr=[hi_radial_fits["peaksub_feath_s_kappa_low_lim"],
                     hi_radial_fits["peaksub_feath_s_kappa_up_lim"]],
               linestyle='-.', drawstyle='steps-mid',
               color=sb.color_palette()[0])

ax[1].set_xlim([0, 8])
ax[1].set_ylim([-0.09, -0.015])
ax[1].grid()
ax[1].set_ylabel(r"$\kappa$")
ax[1].fill_between([0, 0.5], -0.10, -0.01, facecolor='gray', alpha=0.5,
                   zorder=-1)

# Plot the asymm parameters as function of radius
# Horizontal lines for the asymm of the total profiles.
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_asymm'],
               yerr=[hi_radial_fits["peaksub_feath_asymm_low_lim"],
                     hi_radial_fits["peaksub_feath_asymm_up_lim"]],
               linestyle='-', drawstyle='steps-mid',
               color='k', alpha=0.7, label='Total', zorder=-1)
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_n_asymm'],
               yerr=[hi_radial_fits["peaksub_feath_n_asymm_low_lim"],
                     hi_radial_fits["peaksub_feath_n_asymm_up_lim"]],
               linestyle='--', drawstyle='steps-mid',
               color=sb.color_palette()[2], label='North')
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_s_asymm'],
               yerr=[hi_radial_fits["peaksub_feath_s_asymm_low_lim"],
                     hi_radial_fits["peaksub_feath_s_asymm_up_lim"]],
               linestyle='-.', drawstyle='steps-mid',
               color=sb.color_palette()[0], label='South')
# Bold line for the total profile asymmetries
# ax[0].fill_between(hi_radial_fits['bin_center'],
#                    hi_fit_hwhm_vals['Feath. Peak Sub. Params'][-2] -
#                    hi_fit_hwhm_vals['Feath. Peak Sub. Lower Limit'][-2],
#                    hi_fit_hwhm_vals['Feath. Peak Sub. Params'][-2] +
#                    hi_fit_hwhm_vals['Feath. Peak Sub. Upper Limit'][-2],
#                    color='k', #sb.color_palette()[1],
#                    alpha=0.5,
#                    zorder=-1)
# ax[0].fill_between(hi_radial_fits['bin_center'],
#                    hi_fit_hwhm_vals['Feath. Peak Sub. N Params'][-2] -
#                    hi_fit_hwhm_vals['Feath. Peak Sub. N Lower Limit'][-2],
#                    hi_fit_hwhm_vals['Feath. Peak Sub. N Params'][-2] +
#                    hi_fit_hwhm_vals['Feath. Peak Sub. N Upper Limit'][-2],
#                    color=sb.color_palette()[0],
#                    alpha=0.5,
#                    zorder=-1)
# ax[0].fill_between(hi_radial_fits['bin_center'],
#                    hi_fit_hwhm_vals['Feath. Peak Sub. S Params'][-2] -
#                    hi_fit_hwhm_vals['Feath. Peak Sub. S Lower Limit'][-2],
#                    hi_fit_hwhm_vals['Feath. Peak Sub. S Params'][-2] +
#                    hi_fit_hwhm_vals['Feath. Peak Sub. S Upper Limit'][-2],
#                    color=sb.color_palette()[2],
#                    alpha=0.5,
#                    zorder=-1)

ax[0].set_xlim([0, 8])
ax[0].set_ylim([-0.33, 0.22])
ax[0].grid()
ax[0].set_ylabel("Asymmetry")
ax[0].legend(frameon=True)
ax[0].fill_between([0, 0.5], -0.35, 0.25, facecolor='gray', alpha=0.5,
                   zorder=-1)

p.tight_layout()

fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_NS_asymm_kappa_fwing_comparison.pdf")))
fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_NS_asymm_kappa_fwing_comparison.png")))
p.close()

# Plot comparisons of these fits

twocolumn_figure()

fig, ax = p.subplots(2, 3, sharex=True)

bin_cents = hi_radial_fits['bin_center']
hi_params = hi_radial_fits

cpal = sb.color_palette()

ax[0, 0].errorbar(bin_cents, hi_params["peaksub_feath_sigma"],
                  yerr=[hi_params["peaksub_feath_sigma_low_lim"],
                        hi_params["peaksub_feath_sigma_up_lim"]],
                  color='k', label='Total',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[0, 0].errorbar(bin_cents, hi_params["peaksub_feath_n_sigma"],
                  yerr=[hi_params["peaksub_feath_n_sigma_low_lim"],
                        hi_params["peaksub_feath_n_sigma_up_lim"]],
                  color=cpal[0], label='North',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[0, 0].errorbar(bin_cents, hi_params["peaksub_feath_s_sigma"],
                  yerr=[hi_params["peaksub_feath_s_sigma_low_lim"],
                        hi_params["peaksub_feath_s_sigma_up_lim"]],
                  color=cpal[2], label='South',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[0, 0].fill_between([0, 0.5], 5, 9, facecolor='gray', alpha=0.5)
ax[0, 0].set_xlim([0.0, 8.0])
ax[0, 0].grid()
ax[0, 0].set_ylabel(r"$\sigma$ (km/s)")
ax[0, 0].set_ylim([5, 9])
ax[0, 0].legend(frameon=True, loc='upper right')

ax[0, 1].errorbar(bin_cents, hi_params["peaksub_feath_v_peak"],
                  yerr=[hi_params["peaksub_feath_v_peak_low_lim"],
                        hi_params["peaksub_feath_v_peak_up_lim"]],
                  color='k', label='Total',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[0, 1].errorbar(bin_cents, hi_params["peaksub_feath_n_v_peak"],
                  yerr=[hi_params["peaksub_feath_n_v_peak_low_lim"],
                        hi_params["peaksub_feath_n_v_peak_up_lim"]],
                  color=cpal[0], label='North',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[0, 1].errorbar(bin_cents, hi_params["peaksub_feath_s_v_peak"],
                  yerr=[hi_params["peaksub_feath_s_v_peak_low_lim"],
                        hi_params["peaksub_feath_s_v_peak_up_lim"]],
                  color=cpal[2], label='South',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[0, 1].fill_between([0, 0.5], -0.8, 0.8, facecolor='gray', alpha=0.5)
ax[0, 1].set_xlim([0.0, 8.0])
ax[0, 1].grid()
ax[0, 1].set_ylabel(r"$v_{\rm peak}$ (km/s)")
ax[0, 1].set_ylim([-0.8, 0.8])

ax[0, 2].errorbar(bin_cents, hi_params["peaksub_feath_f_wings"],
                  yerr=[hi_params["peaksub_feath_f_wings_low_lim"],
                        hi_params["peaksub_feath_f_wings_up_lim"]],
                  color='k', label='Total',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[0, 2].errorbar(bin_cents, hi_params["peaksub_feath_n_f_wings"],
                  yerr=[hi_params["peaksub_feath_n_f_wings_low_lim"],
                        hi_params["peaksub_feath_n_f_wings_up_lim"]],
                  color=cpal[0], label='North',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[0, 2].errorbar(bin_cents, hi_params["peaksub_feath_s_f_wings"],
                  yerr=[hi_params["peaksub_feath_s_f_wings_low_lim"],
                        hi_params["peaksub_feath_s_f_wings_up_lim"]],
                  color=cpal[2], label='South',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[0, 2].fill_between([0, 0.5], 0.12, 0.37, facecolor='gray', alpha=0.5)
ax[0, 2].set_xlim([0.0, 8.0])
ax[0, 2].grid()
ax[0, 2].set_ylabel(r"$f_{\rm wings}$")
ax[0, 2].set_ylim([0.12, 0.37])

ax[1, 0].errorbar(bin_cents, hi_params["peaksub_feath_sigma_wing"],
                  yerr=[hi_params["peaksub_feath_sigma_wing_low_lim"],
                        hi_params["peaksub_feath_sigma_wing_up_lim"]],
                  color='k', label='Total',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[1, 0].errorbar(bin_cents, hi_params["peaksub_feath_n_sigma_wing"],
                  yerr=[hi_params["peaksub_feath_n_sigma_wing_low_lim"],
                        hi_params["peaksub_feath_n_sigma_wing_up_lim"]],
                  color=cpal[0], label='North',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[1, 0].errorbar(bin_cents, hi_params["peaksub_feath_s_sigma_wing"],
                  yerr=[hi_params["peaksub_feath_s_sigma_wing_low_lim"],
                        hi_params["peaksub_feath_s_sigma_wing_up_lim"]],
                  color=cpal[2], label='South',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[1, 0].fill_between([0, 0.5], 15, 30, facecolor='gray', alpha=0.5)
ax[1, 0].set_xlim([0.0, 8.0])
ax[1, 0].grid()
ax[1, 0].set_xlabel(r"Radius (kpc)")
ax[1, 0].set_ylabel(r"$\sigma_{\rm wings}$ (km/s)")
ax[1, 0].set_ylim([18, 30])

ax[1, 1].errorbar(bin_cents, hi_params["peaksub_feath_asymm"],
                  yerr=[hi_params["peaksub_feath_asymm_low_lim"],
                        hi_params["peaksub_feath_asymm_up_lim"]],
                  color='k', label='Total',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[1, 1].errorbar(bin_cents, hi_params["peaksub_feath_n_asymm"],
                  yerr=[hi_params["peaksub_feath_n_asymm_low_lim"],
                        hi_params["peaksub_feath_n_asymm_up_lim"]],
                  color=cpal[0], label='North',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[1, 1].errorbar(bin_cents, hi_params["peaksub_feath_s_asymm"],
                  yerr=[hi_params["peaksub_feath_s_asymm_low_lim"],
                        hi_params["peaksub_feath_s_asymm_up_lim"]],
                  color=cpal[2], label='South',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[1, 1].fill_between([0, 0.5], -0.22, 0.32, facecolor='gray', alpha=0.5)
ax[1, 1].set_xlim([0.0, 8.0])
ax[1, 1].grid()
ax[1, 1].set_xlabel(r"Radius (kpc)")
ax[1, 1].set_ylabel(r"Asymm.")
ax[1, 1].set_ylim([-0.22, 0.32])

ax[1, 2].errorbar(bin_cents, hi_params["peaksub_feath_kappa"],
                  yerr=[hi_params["peaksub_feath_kappa_low_lim"],
                        hi_params["peaksub_feath_kappa_up_lim"]],
                  color='k', label='Total',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[1, 2].errorbar(bin_cents, hi_params["peaksub_feath_n_kappa"],
                  yerr=[hi_params["peaksub_feath_n_kappa_low_lim"],
                        hi_params["peaksub_feath_n_kappa_up_lim"]],
                  color=cpal[0], label='North',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[1, 2].errorbar(bin_cents, hi_params["peaksub_feath_s_kappa"],
                  yerr=[hi_params["peaksub_feath_s_kappa_low_lim"],
                        hi_params["peaksub_feath_s_kappa_up_lim"]],
                  color=cpal[2], label='South',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[1, 2].fill_between([0, 0.5], -0.09, -0.02, facecolor='gray', alpha=0.5)
ax[1, 2].set_xlim([0.0, 8.0])
ax[1, 2].grid()
ax[1, 2].set_ylabel(r"$\kappa$")
ax[1, 2].set_xlabel(r"Radius (kpc)")
ax[1, 2].set_ylim([-0.09, -0.02])

p.tight_layout()

fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_NS_asymm_comparison_allparams.pdf")))
fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_NS_asymm_comparison_allparams.png")))

p.close()

default_figure()

# There are some interesting differences between the halves. Overplot individual
# stacked profiles

# spec_axis = peakvel_stack_feath.spectral_axis.to(u.km / u.s).value

# for i in range(len(bin_cents)):

#     p.subplot(121)
#     p.plot(spec_axis,
#            peakvel_stack_feath_n[:, i, 0].value / peakvel_stack_feath_n[:, i, 0].max().value,
#            color=cpal[0], drawstyle='steps-mid')
#     p.plot(spec_axis,
#            peakvel_stack_feath_s[:, i, 0].value / peakvel_stack_feath_s[:, i, 0].max().value,
#            color=cpal[2], drawstyle='steps-mid')
#     p.plot(spec_axis,
#            peakvel_stack_feath[:, i, 0].value / peakvel_stack_feath[:, i, 0].max().value,
#            color='k', drawstyle='steps-mid')
#     p.xlim([-50, 50])

#     p.subplot(122)
#     p.plot(spec_axis,
#            peakvel_stack_feath_n[:, i, 0].value / peakvel_stack_feath_n[:, i, 0].max().value -
#            peakvel_stack_feath[:, i, 0].value / peakvel_stack_feath[:, i, 0].max().value,
#            color=cpal[0], drawstyle='steps-mid')
#     p.plot(spec_axis,
#            peakvel_stack_feath_s[:, i, 0].value / peakvel_stack_feath_s[:, i, 0].max().value -
#            peakvel_stack_feath[:, i, 0].value / peakvel_stack_feath[:, i, 0].max().value,
#            color=cpal[2], drawstyle='steps-mid')
#     p.xlim([-50, 50])

#     p.draw()
#     raw_input("{}".format(bin_cents[i]))
#     p.clf()


# Loop over the radial profiles to find the asymmetric line wing fraction

hi_fwing_vals = {}

vels = peakvel_stack_feath_n[:, 0, 0].spectral_axis.to(u.km / u.s).value

for i, b in enumerate(bin_cents):

    # VLA

    params, low_err, up_err = \
        find_linewing_asymm(vels, peakvel_stack_feath_n[:, i, 0].value,
                            peakvel_stack_feath_s[:, i, 0].value,
                            niters=100,
                            sigma_noise_n=sigma_noise * np.sqrt(n_numpix[i] / npix_beam),
                            sigma_noise_s=sigma_noise * np.sqrt(s_numpix[i] / npix_beam))

    hi_fwing_vals[round(b, 2)] = np.hstack([params, low_err, up_err])

hi_fwing_df = DataFrame(hi_fwing_vals).T

onecolumn_figure()

p.errorbar(bin_cents, hi_fwing_df[0],
           yerr=[hi_fwing_df[3], hi_fwing_df[6]],
           label='Total', drawstyle='steps-mid')
p.errorbar(bin_cents, hi_fwing_df[1],
           yerr=[hi_fwing_df[4], hi_fwing_df[7]],
           linestyle='--', label='Symmetric', drawstyle='steps-mid')
p.errorbar(bin_cents, hi_fwing_df[2],
           yerr=[hi_fwing_df[5], hi_fwing_df[8]],
           linestyle='-.', label='Asymmetric', drawstyle='steps-mid')
p.legend(frameon=True, loc='lower right')
p.fill_between([0, 0.5], -0.2, 0.4, facecolor='gray', alpha=0.5,
               zorder=-1)
p.ylim([-0.15, 0.35])
p.xlim([0, 8])
p.grid()
p.ylabel(r"$f_{\rm wings}$")
p.xlabel("Radius (kpc)")
p.tight_layout()

p.savefig(allfigs_path(osjoin(figure_folder, "fwings_NS_symm_asymm.pdf")))
p.savefig(allfigs_path(osjoin(figure_folder, "fwings_NS_symm_asymm.png")))

p.close()

params_tot, low_err_tot, up_err_tot = \
    find_linewing_asymm(vels, total_peakvel_feath_n.value,
                        total_peakvel_feath_s.value,
                        niters=100,
                        sigma_noise_n=sigma_noise * np.sqrt(n_numpix.sum() / npix_beam),
                        sigma_noise_s=sigma_noise * np.sqrt(s_numpix.sum() / npix_beam))

print("Total: {0}+{1}-{2}".format(params_tot[0], low_err_tot[0],
                                  up_err_tot[0]))
print("Symm: {0}+{1}-{2}".format(params_tot[1], low_err_tot[1],
                                 up_err_tot[1]))
print("Asymm: {0}+{1}-{2}".format(params_tot[2], low_err_tot[2],
                                  up_err_tot[2]))
# Total: 0.254678938228+0.00690208853445-0.0100244686344
# Symm: 0.178331243797+0.00702839040532-0.00968523782558
# Asymm: 0.0763476944305+6.35245153629e-05-0.00027544950285

vels = total_peakvel_n.spectral_axis.to(u.km / u.s).value

params_tot_vla, low_err_tot_vla, up_err_tot_vla = \
    find_linewing_asymm(vels, total_peakvel_n.value,
                        total_peakvel_s.value,
                        niters=100,
                        sigma_noise_n=sigma_noise * np.sqrt(n_numpix.sum() / npix_beam),
                        sigma_noise_s=sigma_noise * np.sqrt(s_numpix.sum() / npix_beam))

print("Total: {0}+{1}-{2}".format(params_tot_vla[0], low_err_tot_vla[0],
                                  up_err_tot_vla[0]))
print("Symm: {0}+{1}-{2}".format(params_tot_vla[1], low_err_tot_vla[1],
                                 up_err_tot_vla[1]))
print("Asymm: {0}+{1}-{2}".format(params_tot_vla[2], low_err_tot_vla[2],
                                  up_err_tot_vla[2]))
# Total: 0.186905641003+0.00632970454644-0.0127007065516
# Symm: 0.148370009855+0.00621136896641-0.0122276069357
# Asymm: 0.0385356311475+0.000252790357856-0.000452460576457
