import astropy.units as u
from spectral_cube import OneDSpectrum
import numpy as np
import matplotlib.pyplot as p
from pandas import DataFrame
from astropy.io import fits
from astropy.table import Table
from os.path import join as osjoin
import os

from cube_analysis.spectral_stacking_models import fit_hwhm, fit_gaussian

from paths import (iram_co21_14B088_data_path,
                   allfigs_path, alltables_path, fourteenB_HI_data_wGBT_path,
                   fourteenB_HI_data_path)
# from constants import co21_mass_conversion
from plotting_styles import *


'''
Model and compare the stacked profiles from total_stacked_profiles.py
'''

figure_folder = allfigs_path("stacked_profiles")
if not os.path.exists(figure_folder):
    os.mkdir(figure_folder)


dr = 500 * u.pc
max_radius = (7.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)
maxrad_string = "{0}{1}".format(int(max_radius.value), max_radius.unit)

# Load the CO stacks
co_stackpath = lambda x: osjoin(iram_co21_14B088_data_path("stacked_spectra", no_check=True), x)

total_spectrum_co = OneDSpectrum.from_hdu(fits.open(co_stackpath("rotation_stacked_{}.fits".format(maxrad_string))))
total_spectrum_co_cent = OneDSpectrum.from_hdu(fits.open(co_stackpath("centroid_stacked_{}.fits".format(maxrad_string))))
total_spectrum_co_peakvel = OneDSpectrum.from_hdu(fits.open(co_stackpath("peakvel_stacked_{}.fits".format(maxrad_string))))

# Load the total HI profiles in
comb_stackpath = lambda x: osjoin(fourteenB_HI_data_wGBT_path("", no_check=True), "stacked_spectra", x)

total_spectrum_hi = OneDSpectrum.from_hdu(fits.open(comb_stackpath("rotation_stacked.fits")))
total_spectrum_hi_cent = OneDSpectrum.from_hdu(fits.open(comb_stackpath("centroid_stacked.fits")))
total_spectrum_hi_peakvel = OneDSpectrum.from_hdu(fits.open(comb_stackpath("peakvel_stacked.fits")))

# Define the shifted CO velocity axis
twocolumn_twopanel_figure()
# Plot the profiles.
fig, ax = p.subplots(1, 3, sharey=True, sharex=True)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)
ax[0].plot(total_spectrum_hi.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi / total_spectrum_hi.max()).value,
           '-', drawstyle='steps-mid', label="HI")
ax[0].plot(total_spectrum_co.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co / total_spectrum_co.max()).value,
           '--', drawstyle='steps-mid', label="CO(2-1)")
ax[0].set_xlabel("Velocity (km/s)")
ax[0].set_ylabel("Normalized Total Intensity")
# ax[0].set_title("Rotation subtracted")
ax[0].text(-55, 0.88, "Rotation\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].grid()

ax[1].plot(total_spectrum_hi_cent.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi_cent / total_spectrum_hi_cent.max()).value,
           '-', drawstyle='steps-mid')
ax[1].plot(total_spectrum_co_cent.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co_cent / total_spectrum_co_cent.max()).value,
           '--', drawstyle='steps-mid')
# ax[1].set_title("Centroid subtracted")
ax[1].text(-55, 0.88, "Centroid\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[1].set_xlabel("Velocity (km/s)")
ax[1].grid()

ax[2].plot(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi_peakvel / total_spectrum_hi_peakvel.max()).value,
           '-', drawstyle='steps-mid', label="HI")
ax[2].plot(total_spectrum_co_peakvel.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co_peakvel / total_spectrum_co_peakvel.max()).value,
           '--', drawstyle='steps-mid', label="CO(2-1)")
# ax[2].set_title("Centroid subtracted")
ax[2].text(-55, 0.88, "Peak Vel.\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].set_xlabel("Velocity (km/s)")
ax[2].set_ylim([-0.02, 1.1])
ax[2].set_xlim([-60, 60])
ax[2].legend(frameon=True, fontsize=9)
ax[2].grid()

p.tight_layout()

fig.savefig(osjoin(figure_folder, "total_profile_corrected_velocity_HI_CO21.pdf"))
fig.savefig(osjoin(figure_folder, "total_profile_corrected_velocity_HI_CO21.png"))

p.close()
# raw_input("Next plot?")

# Total CO mass. Using 6.7 Msol / pc^2 / K km s^-1\
# pixscale = gal.distance.to(u.pc) * (np.pi / 180.) * \
#     np.abs(co_cube.header["CDELT2"])
# chan_width = \
#     np.abs(co_cube.spectral_axis[1] - co_cube.spectral_axis[0]).to(u.km / u.s)
# beam_eff = 0.75  # Beam efficiency of IRAM @ 235 GHz
# # Where total_spectrum_co is in K
# total_co_mass = \
#     total_spectrum_co[total_spectrum_co > 0 * u.K].sum() * chan_width * \
#     pixscale ** 2 * co21_mass_conversion / beam_eff
# print("Total H2 Mass from CO is {} Msol".format(total_co_mass))


# Model the width, based on scaling from the HWHM

# Load in HI fits for the total profiles
hi_tab = Table.read(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_hwhm_fits.csv"))

# CO model

spectra = [total_spectrum_co, total_spectrum_co_cent,
           total_spectrum_co_peakvel]
co_fit_vals = {}
hwhm_models = {}
labels = ["Rot. Sub.", "Cent. Sub.", "Peak Sub."]
file_labels = ["rotsub", "centsub", "peaksub"]

sigma_noise = (20.33 * u.mK).to(u.K).value
num_pix_total = (fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_source_mask.fits"))[0].data > 0).sum()

# Data is smoothed and reprojected to match the HI beam size
npix_beam = 41.

onecolumn_figure()

for spectrum, label, file_label in zip(spectra, labels, file_labels):

    norm_intens = spectrum / spectrum.max()
    vels = spectrum.spectral_axis.to(u.km / u.s).value

    # HWHM fitting
    parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_CO = \
        fit_hwhm(vels, norm_intens,
                 sigma_noise=sigma_noise / np.nanmax(spectrum.value),
                 nbeams=num_pix_total / npix_beam)

    # parvals_gauss = fit_gaussian(vels, norm_intens)
    # print(parvals_hwhm[0], parvals_gauss[0])

    co_fit_vals[label + " Params"] = parvals_hwhm
    co_fit_vals[label + " Errors"] = parerrs_hwhm
    hwhm_models[file_label] = g_CO

    # Better sampling for plotting
    more_vels = np.arange(vels.min(), vels.max(), 0.5)

    ax1 = p.subplot2grid((5, 1), (0, 0), rowspan=4)
    p.plot(vels, norm_intens, 'b-', drawstyle='steps-mid')
    p.plot(more_vels, g_CO(more_vels), 'k--', label="Fit")
    p.ylabel("CO Normalized Intensity")
    p.xticks([])
    p.xlim([-50, 50])
    p.legend(frameon=True)
    p.ylim([-0.1, 1.1])
    p.grid()

    ax2 = p.subplot2grid((5, 1), (4, 0))
    p.plot(vels, norm_intens - g_CO(vels), 'b-',
           drawstyle='steps-mid')
    p.hlines(0.0, -50, 50, color='k')
    p.grid()
    # p.ylim([-0.015, 0.015])
    p.xlim([-50, 50])
    p.yticks([-0.06, -0.03, 0.0, 0.03, 0.06])
    p.xlabel("Velocity (km/s)")
    p.ylabel("Residuals")

    p.tight_layout()
    p.subplots_adjust(hspace=0)

    filename = "total_profile_stack_{}_co21_hwhm_fit".format(file_label)
    p.savefig(osjoin(figure_folder, filename + ".pdf"))
    p.savefig(osjoin(figure_folder, filename + ".png"))
    # p.draw()
    # raw_input("Next plot?")
    # p.clf()
    p.close()

# Save table of parameters
co_param_df = DataFrame(co_fit_vals, index=parnames_hwhm)
co_param_df.to_latex(alltables_path("co_gaussian_totalprof_fits_hwhm.tex"))
co_param_df.to_csv(iram_co21_14B088_data_path("tables/co_gaussian_totalprof_fits_hwhm.csv",
                                              no_check=True))

# Load in radial profiles

total_spectrum_hi_radial = SpectralCube.read(comb_stackpath("rotation_stacked_radial_{}.fits".format(wstring)))
total_spectrum_hi_radial_cent = SpectralCube.read(comb_stackpath("centroid_stacked_radial_{}.fits".format(wstring)))
total_spectrum_hi_radial_peakvel = SpectralCube.read(comb_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)))

total_spectrum_co_radial = SpectralCube.read(co_stackpath("rotation_stacked_radial_{}.fits".format(wstring)))
total_spectrum_co_radial_cent = SpectralCube.read(co_stackpath("centroid_stacked_radial_{}.fits".format(wstring)))
total_spectrum_co_radial_peakvel = SpectralCube.read(co_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)))

# Per radial bin spectra
Nrows = 4
Ncols = 4

twocolumn_figure(font_scale=1.3)

co_radstacks = [total_spectrum_co_radial, total_spectrum_co_radial_cent,
                total_spectrum_co_radial_peakvel]
hi_radstacks = [total_spectrum_hi_radial, total_spectrum_hi_radial_cent,
                total_spectrum_hi_radial_peakvel]
labels = ["rotsub", "centsub", "peaksub"]

# Make the bin edges
# Convert these into kpc
inneredge = np.arange(0, max_radius.value, dr.value) / 1000.
outeredge = (inneredge + dr.to(u.kpc).value)

for co_stack, hi_stack, label in zip(co_radstacks, hi_radstacks, labels):
    p.figure(1, figsize=(8.4, 11)).clf()

    fig, ax = p.subplots(Nrows, Ncols,
                         sharex=True,
                         sharey=True, num=1)

    p.subplots_adjust(hspace=0.05,
                      wspace=0.05)

    fig.text(0.5, 0.04, 'Velocity (km/s)', ha='center')
    fig.text(0.04, 0.5, 'Normalized Intensity', va='center',
             rotation='vertical')

    for ctr, (r0, r1) in enumerate(zip(inneredge,
                                       outeredge)):

        r, c = np.unravel_index(ctr, (Nrows, Ncols))

        norm_hi = (hi_stack[:, ctr, 0] /
                   hi_stack[:, ctr, 0].max()).value

        norm_co = (co_stack[:, ctr, 0] /
                   co_stack[:, ctr, 0].max()).value

        ax[r, c].plot(hi_stack.spectral_axis.to(u.km / u.s).value,
                      norm_hi,
                      '-', drawstyle='steps-mid', label="HI", alpha=0.7)
        # There's a 1 channel offset from my rotation subtraction in the cube
        ax[r, c].plot(co_stack.spectral_axis.to(u.km / u.s).value,
                      norm_co,
                      '--', drawstyle='steps-mid', label="CO(2-1)", alpha=0.7)
        ax[r, c].set_ylim([-0.02, 1.1])
        ax[r, c].set_xlim([-70, 70])

        ax[r, c].annotate("{0}-{1} kpc".format(r0, r1),
                          xy=(13, 0.95),
                          color='k', fontsize=8,
                          bbox={"boxstyle": "square", "facecolor": "w"})

        if ctr == 0:
            ax[r, c].legend(loc='lower left', frameon=True, prop={"size": 8})
        ax[r, c].grid()

    fig.savefig(osjoin(figure_folder, "total_profile_{}_hi_co_radial.pdf".format(label)))
    # fig.savefig(osjoin(figure_folder, "total_profile_{}_hi_co_radial.png".format(label)))

    p.close()

# How do the model parameters change with radius?

# Pixels in each radial bin
num_pix = np.load(iram_co21_14B088_data_path("stacked_spectra/radial_stacking_pixelsinbin_{}.npy".format(wstring), no_check=True))

hi_params = {}
co_params = {}

param_names = ["sigma", "v_peak", "f_wings", "sigma_wing", "asymm", "kappa"]

sigma_noise_hi = 2.8  # K

for sub in labels:
    for name in param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_error = "{}_stderr".format(par_name)

        hi_params[par_name] = np.zeros_like(inneredge)
        hi_params[par_error] = np.zeros_like(inneredge)

        co_params[par_name] = np.zeros_like(inneredge)
        co_params[par_error] = np.zeros_like(inneredge)


for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    hi_spectra = [total_spectrum_hi_radial[:, ctr, 0],
                  total_spectrum_hi_radial_cent[:, ctr, 0],
                  total_spectrum_hi_radial_peakvel[:, ctr, 0]]

    for spectrum, label in zip(hi_spectra, labels):

        vels = spectrum.spectral_axis.to(u.km / u.s).value
        norm_intens = (spectrum / spectrum.max()).value

        params, stderrs, names = \
            fit_hwhm(vels, spectrum, sigma_noise=sigma_noise_hi,
                     nbeams=num_pix[ctr] / npix_beam)[:-1]

        for idx, name in enumerate(names):
            par_name = "{0}_{1}".format(label, name)
            hi_params[par_name][ctr] = params[idx]
            hi_params["{}_stderr".format(par_name)][ctr] = stderrs[idx]

    co_spectra = [total_spectrum_co_radial[:, ctr, 0],
                  total_spectrum_co_radial_cent[:, ctr, 0],
                  total_spectrum_co_radial_peakvel[:, ctr, 0]]

    for spectrum, label in zip(co_spectra, labels):

        vels = spectrum.spectral_axis.to(u.km / u.s).value
        norm_intens = (spectrum / spectrum.max()).value

        params, stderrs, names = \
            fit_hwhm(vels, spectrum, sigma_noise=sigma_noise,
                     nbeams=num_pix[ctr] / npix_beam)[:-1]

        for idx, name in enumerate(names):
            par_name = "{0}_{1}".format(label, name)
            co_params[par_name][ctr] = params[idx]
            co_params["{}_stderr".format(par_name)][ctr] = stderrs[idx]

bin_names = ["{}-{}".format(r0, r1)
             for r0, r1 in zip(inneredge, outeredge)]

co_radial_fits = DataFrame(co_params, index=bin_names)
hi_radial_fits = DataFrame(hi_params, index=bin_names)

co_radial_fits.to_latex(alltables_path("co_hwhm_totalprof_fits_radial_{}.tex".format(wstring)))
co_radial_fits.to_csv(iram_co21_14B088_data_path("tables/co_hwhm_totalprof_fits_radial_{}.csv".format(wstring),
                                                 no_check=True))

hi_radial_fits.to_latex(alltables_path("hi_hwhm_totalprof_fits_radial_{}.tex".format(wstring)))
hi_radial_fits.to_csv(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_radial_{}.csv".format(wstring),
                                             no_check=True))

# Plot comparisons of these fits
bin_cents = outeredge - dr.to(u.kpc).value / 2.

twocolumn_twopanel_figure()

fig, ax = p.subplots(1, 3, sharey=True)

ax[0].errorbar(bin_cents, hi_params["rotsub_sigma"],
               yerr=hi_params["rotsub_sigma_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_sigma"],
               yerr=co_params["rotsub_sigma_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].legend(loc='lower left', frameon=True)
ax[0].grid()
ax[0].set_ylim([0.25, 16])
ax[0].text(1.3, 14, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel("Gaussian Width (km/s)")

ax[1].errorbar(bin_cents, hi_params["centsub_sigma"],
               yerr=hi_params["centsub_sigma_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_sigma"],
               yerr=co_params["centsub_sigma_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(1.3, 14, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].errorbar(bin_cents, hi_params["peaksub_sigma"],
               yerr=hi_params["peaksub_sigma_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_sigma"],
               yerr=co_params["peaksub_sigma_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
ax[2].set_xlabel("Radius (kpc)")
ax[2].text(1.3, 14, "Peak Vel. subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
p.tight_layout()
p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.savefig(osjoin(figure_folder, "total_profile_radial_widths_HI_CO21.pdf"))
fig.savefig(osjoin(figure_folder, "total_profile_radial_widths_HI_CO21.png"))

p.close()

# Compare velocity at the peak of the stacks

fig, ax = p.subplots(1, 3, sharey=True)

ax[0].errorbar(bin_cents, hi_params["rotsub_v_peak"],
               yerr=hi_params["rotsub_v_peak_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_v_peak"],
               yerr=co_params["rotsub_v_peak_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].legend(loc='lower left', frameon=True)
ax[0].grid()
# ax[0].set_ylim([0.25, 16])
ax[0].text(1.3, 10, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel("Central Velocity (km/s)")

ax[1].errorbar(bin_cents, hi_params["centsub_v_peak"],
               yerr=hi_params["centsub_v_peak_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_v_peak"],
               yerr=co_params["centsub_v_peak_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(1.3, 10, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].errorbar(bin_cents, hi_params["peaksub_v_peak"],
               yerr=hi_params["peaksub_v_peak_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_v_peak"],
               yerr=co_params["peaksub_v_peak_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
ax[2].set_xlabel("Radius (kpc)")
ax[2].text(1.3, 10, "Peak Vel. subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
p.tight_layout()
p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.savefig(osjoin(figure_folder, "total_profile_radial_vpeak_HI_CO21.pdf"))
fig.savefig(osjoin(figure_folder, "total_profile_radial_vpeak_HI_CO21.png"))

p.close()

# Compare f_wings

fig, ax = p.subplots(1, 3, sharey=True)

ax[0].errorbar(bin_cents, hi_params["rotsub_f_wings"],
               yerr=hi_params["rotsub_f_wings_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_f_wings"],
               yerr=co_params["rotsub_f_wings_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].legend(loc='lower left', frameon=True)
ax[0].grid()
ax[0].set_ylim([-1, 1])
ax[0].text(1.3, 0.8, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel(r"$f_{\rm wings}$")

ax[1].errorbar(bin_cents, hi_params["centsub_f_wings"],
               yerr=hi_params["centsub_f_wings_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_f_wings"],
               yerr=co_params["centsub_f_wings_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(1.3, 0.8, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].errorbar(bin_cents, hi_params["peaksub_f_wings"],
               yerr=hi_params["peaksub_f_wings_stderr"],
               label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_f_wings"],
               yerr=co_params["peaksub_f_wings_stderr"],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
ax[2].set_xlabel("Radius (kpc)")
ax[2].text(1.3, 0.8, "Peak Vel. subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
p.tight_layout()
p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.savefig(osjoin(figure_folder, "total_profile_radial_fwings_HI_CO21.pdf"))
fig.savefig(osjoin(figure_folder, "total_profile_radial_fwings_HI_CO21.png"))

p.close()

# The other parameters don't have interesting features. sigma_wing is highly
# uncertain. The asymmetry is affected by the large velocity channels and I
# don't trust it. kappa are all consistent with a Gaussian shape.

p.close()
default_figure()