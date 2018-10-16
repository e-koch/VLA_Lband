import astropy.units as u
import astropy.constants as const
from spectral_cube import OneDSpectrum, SpectralCube
import numpy as np
import matplotlib.pyplot as p
from pandas import DataFrame
from astropy.io import fits
from astropy.table import Table
from os.path import join as osjoin
import os
import seaborn as sb
from astropy.modeling.models import Gaussian1D

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

# fig.savefig(osjoin(figure_folder, "total_profile_corrected_velocity_HI_CO21.pdf"))
# fig.savefig(osjoin(figure_folder, "total_profile_corrected_velocity_HI_CO21.png"))

p.close()

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

sigma_noise = (16 * u.mK).to(u.K).value
num_pix_total = (fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_source_mask.fits"))[0].data > 0).sum()

# Data is smoothed and reprojected to match the HI beam size
npix_beam = 41.

onecolumn_figure()

for spectrum, label, file_label in zip(spectra, labels, file_labels):

    norm_intens = spectrum / spectrum.max()
    vels = spectrum.spectral_axis.to(u.km / u.s).value

    # Restrict to +/- 50 km/s.
    # CO data has some nasty systematics at large radii
    vel_mask = np.logical_and(vels < 50, vels > -50)

    # HWHM fitting
    parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_CO = \
        fit_hwhm(vels[vel_mask], norm_intens[vel_mask],
                 sigma_noise=sigma_noise / np.nanmax(spectrum.value),
                 nbeams=num_pix_total / npix_beam,
                 niters=1000)[:-1]

    # parvals_gauss = fit_gaussian(vels, norm_intens)
    # print(parvals_hwhm[0], parvals_gauss[0])

    co_fit_vals[label + " Params"] = parvals_hwhm
    co_fit_vals[label + " Lower Limit"] = parerrs_hwhm[0]
    co_fit_vals[label + " Upper Limit"] = parerrs_hwhm[1]
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

# Load up already saved version
# import pandas as pd
# co_param_df = pd.read_csv(iram_co21_14B088_data_path("tables/co_gaussian_totalprof_fits_hwhm.csv"), index_col=0)

# Single plot of just the peak velocity profiles with the width highlighted
onecolumn_figure()
cpal = sb.color_palette()
p.plot(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value,
       (total_spectrum_hi_peakvel / total_spectrum_hi_peakvel.max()).value,
       '-', drawstyle='steps-mid', label="HI", color=cpal[2])
p.plot(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value,
       Gaussian1D(1.0, hi_tab['Feath. Peak Sub. Params'][1],
                  hi_tab['Feath. Peak Sub. Params'][0])(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value),
       color=cpal[2], linewidth=4, alpha=0.5, zorder=-1)
# p.axvline(hi_tab['Feath. Peak Sub. Params'][1] +
#           hi_tab['Feath. Peak Sub. Params'][0] *
#           np.sqrt(2 * np.log(2)), linestyle='-', color=cpal[2], linewidth=3,
#           alpha=0.75)
# p.axvline(hi_tab['Feath. Peak Sub. Params'][1] -
#           hi_tab['Feath. Peak Sub. Params'][0] *
#           np.sqrt(2 * np.log(2)), linestyle='-', color=cpal[2], linewidth=3,
#           alpha=0.75)
p.plot(total_spectrum_co_peakvel.spectral_axis.to(u.km / u.s).value,
       (total_spectrum_co_peakvel / total_spectrum_co_peakvel.max()).value,
       '--', drawstyle='steps-mid', label="CO(2-1)", color=cpal[1])
p.plot(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value,
       Gaussian1D(1.0, co_param_df['Peak Sub. Params'][1],
                  co_param_df['Peak Sub. Params'][0])(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value),
       color=cpal[1], linewidth=4, alpha=0.5, zorder=-1)

# p.axvline(co_param_df['Peak Sub. Params'][1] +
#           co_param_df['Peak Sub. Params'][0] *
#           np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1], linewidth=3,
#           alpha=0.75)
# p.axvline(co_param_df['Peak Sub. Params'][1] -
#           co_param_df['Peak Sub. Params'][0] *
#           np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1], linewidth=3,
#           alpha=0.75)
p.ylabel("Normalized Total Intensity")
p.xlabel("Velocity (km/s)")
p.ylim([-0.02, 1.1])
p.xlim([-50, 50])
p.legend(frameon=True, fontsize=9)
p.grid()

p.tight_layout()

p.savefig(osjoin(figure_folder, "total_profile_peak_velocity_HI_CO21.pdf"))
p.savefig(osjoin(figure_folder, "total_profile_peak_velocity_HI_CO21.png"))

p.close()

# And with the equivalent Gaussian models

onecolumn_Npanel_figure(2.5)

vels = total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value

fig, ax = p.subplots(3, 1, sharex=True, sharey=True)

ax[0].plot(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi_peakvel / total_spectrum_hi_peakvel.max()).value,
           '-', drawstyle='steps-mid', label="HI")
ax[0].axvline(hi_tab['Feath. Peak Sub. Params'][1] +
              hi_tab['Feath. Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='-', color=cpal[2],
              linewidth=3,
              alpha=0.75)
ax[0].axvline(hi_tab['Feath. Peak Sub. Params'][1] -
              hi_tab['Feath. Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='-', color=cpal[2],
              linewidth=3,
              alpha=0.75)
ax[0].plot(total_spectrum_co_peakvel.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co_peakvel / total_spectrum_co_peakvel.max()).value,
           '--', drawstyle='steps-mid', label="CO(2-1)")
ax[0].axvline(co_param_df['Peak Sub. Params'][1] +
              co_param_df['Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1],
              linewidth=3,
              alpha=0.75)
ax[0].axvline(co_param_df['Peak Sub. Params'][1] -
              co_param_df['Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1],
              linewidth=3,
              alpha=0.75)
ax[0].set_ylabel("Normalized Total Intensity")
ax[0].grid()

ax[1].plot(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi_peakvel / total_spectrum_hi_peakvel.max()).value,
           '-', drawstyle='steps-mid', label="HI")
ax[1].plot(vels,
           Gaussian1D(1.0, hi_tab['Feath. Peak Sub. Params'][1],
                      hi_tab['Feath. Peak Sub. Params'][0])(vels),
           color='gray',
           alpha=0.5, linewidth=4)
ax[1].axvline(hi_tab['Feath. Peak Sub. Params'][1] +
              hi_tab['Feath. Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='-', color=cpal[2],
              linewidth=3,
              alpha=0.75)
ax[1].axvline(hi_tab['Feath. Peak Sub. Params'][1] -
              hi_tab['Feath. Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='-', color=cpal[2],
              linewidth=3,
              alpha=0.75)
ax[1].set_ylabel("Normalized Total Intensity")
ax[1].grid()
ax[1].text(35, 0.95, "HI", ha='center', va='center',
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].plot(total_spectrum_co_peakvel.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co_peakvel / total_spectrum_co_peakvel.max()).value,
           '--', drawstyle='steps-mid', label="CO(2-1)", color=cpal[1])
ax[2].plot(vels,
           Gaussian1D(1.0, co_param_df['Peak Sub. Params'][1],
                      co_param_df['Peak Sub. Params'][0])(vels),
           color='gray',
           alpha=0.5, linewidth=4)
ax[2].axvline(co_param_df['Peak Sub. Params'][1] +
              co_param_df['Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1],
              linewidth=3,
              alpha=0.75)
ax[2].axvline(co_param_df['Peak Sub. Params'][1] -
              co_param_df['Peak Sub. Params'][0] *
              np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1],
              linewidth=3,
              alpha=0.75)
ax[2].set_ylabel("Normalized Total Intensity")
ax[2].grid()
ax[2].text(35, 0.95, "CO(2-1)", ha='center', va='center',
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].set_xlabel("Velocity (km/s)")
ax[2].set_ylim([-0.02, 1.1])
ax[2].set_xlim([-50, 50])
ax[0].legend(frameon=True, fontsize=9)

p.tight_layout()

p.savefig(osjoin(figure_folder, "total_profile_peak_velocity_w_model_HI_CO21.pdf"))
p.savefig(osjoin(figure_folder, "total_profile_peak_velocity_w_model_HI_CO21.png"))

p.close()

# Load in radial profiles

total_spectrum_hi_radial = SpectralCube.read(comb_stackpath("rotation_stacked_radial_{}.fits".format(wstring)))
total_spectrum_hi_radial_n = SpectralCube.read(comb_stackpath("rotation_stacked_radial_north_{}.fits".format(wstring)))
total_spectrum_hi_radial_s = SpectralCube.read(comb_stackpath("rotation_stacked_radial_south_{}.fits".format(wstring)))

total_spectrum_hi_radial_cent = SpectralCube.read(comb_stackpath("centroid_stacked_radial_{}.fits".format(wstring)))
total_spectrum_hi_radial_cent_n = SpectralCube.read(comb_stackpath("centroid_stacked_radial_north_{}.fits".format(wstring)))
total_spectrum_hi_radial_cent_s = SpectralCube.read(comb_stackpath("centroid_stacked_radial_south_{}.fits".format(wstring)))

total_spectrum_hi_radial_peakvel = SpectralCube.read(comb_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)))
total_spectrum_hi_radial_peakvel_n = SpectralCube.read(comb_stackpath("peakvel_stacked_radial_north_{}.fits".format(wstring)))
total_spectrum_hi_radial_peakvel_s = SpectralCube.read(comb_stackpath("peakvel_stacked_radial_south_{}.fits".format(wstring)))

total_spectrum_co_radial = SpectralCube.read(co_stackpath("rotation_stacked_radial_{}.fits".format(wstring)))
total_spectrum_co_radial_n = SpectralCube.read(co_stackpath("rotation_stacked_radial_north_{}.fits".format(wstring)))
total_spectrum_co_radial_s = SpectralCube.read(co_stackpath("rotation_stacked_radial_south_{}.fits".format(wstring)))

total_spectrum_co_radial_cent = SpectralCube.read(co_stackpath("centroid_stacked_radial_{}.fits".format(wstring)))
total_spectrum_co_radial_cent_n = SpectralCube.read(co_stackpath("centroid_stacked_radial_north_{}.fits".format(wstring)))
total_spectrum_co_radial_cent_s = SpectralCube.read(co_stackpath("centroid_stacked_radial_south_{}.fits".format(wstring)))

total_spectrum_co_radial_peakvel = SpectralCube.read(co_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)))
total_spectrum_co_radial_peakvel_n = SpectralCube.read(co_stackpath("peakvel_stacked_radial_north_{}.fits".format(wstring)))
total_spectrum_co_radial_peakvel_s = SpectralCube.read(co_stackpath("peakvel_stacked_radial_south_{}.fits".format(wstring)))

labels = ["rotsub", "rotsub_n", "rotsub_s",
          "centsub", "centsub_n", "centsub_s",
          "peaksub", "peaksub_n", "peaksub_s"]

# Make the bin edges
# Convert these into kpc
inneredge = np.arange(0, max_radius.value, dr.value) / 1000.
outeredge = (inneredge + dr.to(u.kpc).value)

# How do the model parameters change with radius?

# Pixels in each radial bin
num_pix = np.load(iram_co21_14B088_data_path("stacked_spectra/radial_stacking_pixelsinbin_{}.npy".format(wstring), no_check=True))
num_pix_n = np.load(iram_co21_14B088_data_path("stacked_spectra/radial_stacking_pixelsinbin_north_{}.npy".format(wstring), no_check=True))
num_pix_s = np.load(iram_co21_14B088_data_path("stacked_spectra/radial_stacking_pixelsinbin_south_{}.npy".format(wstring), no_check=True))

hi_params = {}
co_params = {}

hi_models = {}
co_models = {}

param_names = ["sigma", "v_peak", "f_wings", "sigma_wing", "asymm", "kappa"]

sigma_noise_hi = 2.8  # K

for sub in labels:
    for name in param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_lowlim = "{}_low_lim".format(par_name)
        par_uplim = "{}_up_lim".format(par_name)

        hi_params[par_name] = np.zeros_like(inneredge)
        hi_params[par_lowlim] = np.zeros_like(inneredge)
        hi_params[par_uplim] = np.zeros_like(inneredge)

        co_params[par_name] = np.zeros_like(inneredge)
        co_params[par_lowlim] = np.zeros_like(inneredge)
        co_params[par_uplim] = np.zeros_like(inneredge)

    hi_models[sub] = []
    co_models[sub] = []

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):
    print(ctr, len(inneredge))
    hi_spectra = [total_spectrum_hi_radial[:, ctr, 0],
                  total_spectrum_hi_radial_n[:, ctr, 0],
                  total_spectrum_hi_radial_s[:, ctr, 0],
                  total_spectrum_hi_radial_cent[:, ctr, 0],
                  total_spectrum_hi_radial_cent_n[:, ctr, 0],
                  total_spectrum_hi_radial_cent_s[:, ctr, 0],
                  total_spectrum_hi_radial_peakvel[:, ctr, 0],
                  total_spectrum_hi_radial_peakvel_n[:, ctr, 0],
                  total_spectrum_hi_radial_peakvel_s[:, ctr, 0]]

    if "_n" in label:
        pix_in_bin = num_pix_n[ctr]
    elif "_s" in label:
        pix_in_bin = num_pix_s[ctr]
    else:
        pix_in_bin = num_pix[ctr]

    for spectrum, label in zip(hi_spectra, labels):

        vels = spectrum.spectral_axis.to(u.km / u.s).value
        norm_intens = (spectrum / spectrum.max()).value

        vel_mask = np.logical_and(vels < 50, vels > -50)

        params, stderrs, names, hwhm_mod = \
            fit_hwhm(vels[vel_mask], spectrum[vel_mask],
                     sigma_noise=sigma_noise_hi,
                     nbeams=pix_in_bin / npix_beam,
                     niters=1000)[:-1]
        hi_models[label].append(hwhm_mod)

        for idx, name in enumerate(names):
            par_name = "{0}_{1}".format(label, name)
            hi_params[par_name][ctr] = params[idx]
            hi_params["{}_low_lim".format(par_name)][ctr] = \
                np.abs(stderrs[0, idx])
            hi_params["{}_up_lim".format(par_name)][ctr] = \
                np.abs(stderrs[1, idx])

    co_spectra = [total_spectrum_co_radial[:, ctr, 0],
                  total_spectrum_co_radial_n[:, ctr, 0],
                  total_spectrum_co_radial_s[:, ctr, 0],
                  total_spectrum_co_radial_cent[:, ctr, 0],
                  total_spectrum_co_radial_cent_n[:, ctr, 0],
                  total_spectrum_co_radial_cent_s[:, ctr, 0],
                  total_spectrum_co_radial_peakvel[:, ctr, 0],
                  total_spectrum_co_radial_peakvel_n[:, ctr, 0],
                  total_spectrum_co_radial_peakvel_s[:, ctr, 0]]

    for spectrum, label in zip(co_spectra, labels):

        vels = spectrum.spectral_axis.to(u.km / u.s).value
        norm_intens = (spectrum / spectrum.max()).value

        vel_mask = np.logical_and(vels < 50, vels > -50)

        params, stderrs, names, hwhm_mod = \
            fit_hwhm(vels[vel_mask], spectrum[vel_mask],
                     sigma_noise=sigma_noise,
                     nbeams=num_pix[ctr] / npix_beam,
                     niters=1000)[:-1]

        co_models[label].append(hwhm_mod)

        for idx, name in enumerate(names):
            par_name = "{0}_{1}".format(label, name)
            co_params[par_name][ctr] = params[idx]
            co_params["{}_low_lim".format(par_name)][ctr] = \
                np.abs(stderrs[0, idx])
            co_params["{}_up_lim".format(par_name)][ctr] = \
                np.abs(stderrs[1, idx])

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

# Load previously saved version
# co_radial_fits = pd.read_csv(iram_co21_14B088_data_path("tables/co_hwhm_totalprof_fits_radial_{}.csv".format(wstring)),
#                              index_col=0)
# hi_radial_fits = pd.read_csv(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_radial_{}.csv".format(wstring)),
#                              index_col=0)
# co_params = co_radial_fits
# hi_params = hi_radial_fits

# Overplot HI and CO stacks with their equivalent models

co_radstacks = [total_spectrum_co_radial,
                total_spectrum_co_radial_n,
                total_spectrum_co_radial_s,
                total_spectrum_co_radial_cent,
                total_spectrum_co_radial_cent_n,
                total_spectrum_co_radial_cent_s,
                total_spectrum_co_radial_peakvel,
                total_spectrum_co_radial_peakvel_n,
                total_spectrum_co_radial_peakvel_s]

hi_radstacks = [total_spectrum_hi_radial,
                total_spectrum_hi_radial_n,
                total_spectrum_hi_radial_s,
                total_spectrum_hi_radial_cent,
                total_spectrum_hi_radial_cent_n,
                total_spectrum_hi_radial_cent_s,
                total_spectrum_hi_radial_peakvel,
                total_spectrum_hi_radial_peakvel_n,
                total_spectrum_hi_radial_peakvel_s]

# Per radial bin spectra
Nrows = 4
Ncols = 4

twocolumn_figure(font_scale=1.3)

cpal = sb.color_palette()

for co_stack, hi_stack, label in zip(co_radstacks, hi_radstacks, labels):
    p.figure(1, figsize=(8.4, 11)).clf()

    # Grab the HWHM models
    co_mods = co_models[label]
    hi_mods = hi_models[label]

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

        vels = np.linspace(-70, 70, 300)

        ax[r, c].plot(vels, hi_mods[ctr](vels) / hi_stack[:, ctr, 0].max().value,
                      '-', color=cpal[0], alpha=0.5, linewidth=3,)
        ax[r, c].plot(vels, co_mods[ctr](vels) / co_stack[:, ctr, 0].max().value,
                      '-', color=cpal[1], alpha=0.5, linewidth=3,)

        ax[r, c].plot(hi_stack.spectral_axis.to(u.km / u.s).value,
                      norm_hi, '-', color=cpal[0],
                      drawstyle='steps-mid', label="HI", alpha=0.9)
        ax[r, c].plot(co_stack.spectral_axis.to(u.km / u.s).value,
                      norm_co,
                      '--', drawstyle='steps-mid', label="CO(2-1)", alpha=0.9)
        ax[r, c].set_ylim([-0.02, 1.1])
        ax[r, c].set_xlim([-70, 70])

        ax[r, c].annotate("{0}-{1} kpc".format(r0, r1),
                          xy=(13, 0.95),
                          color='k', fontsize=8,
                          bbox={"boxstyle": "square", "facecolor": "w"})

        if ctr == 0:
            ax[r, c].legend(loc='lower left', frameon=True, prop={"size": 8})
        ax[r, c].grid()

    # print(label)
    # raw_input("?")

    fig.savefig(osjoin(figure_folder, "total_profile_{}_hi_co_radial.pdf".format(label)))
    # fig.savefig(osjoin(figure_folder, "total_profile_{}_hi_co_radial.png".format(label)))

    p.close()

# Plot comparisons of these fits
bin_cents = outeredge - dr.to(u.kpc).value / 2.

twocolumn_twopanel_figure()

fig, ax = p.subplots(1, 3, sharey=True)

ax[0].errorbar(bin_cents, hi_params["rotsub_sigma"],
               yerr=[hi_params["rotsub_sigma_low_lim"],
                     hi_params["rotsub_sigma_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_sigma"],
               yerr=[co_params["rotsub_sigma_low_lim"],
                     co_params["rotsub_sigma_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].legend(loc='lower left', frameon=True)
ax[0].grid()
ax[0].set_ylim([0.25, 16])
ax[0].text(1.3, 14, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
# ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel("Gaussian Width (km/s)")

ax[1].errorbar(bin_cents, hi_params["centsub_sigma"],
               yerr=[hi_params["centsub_sigma_low_lim"],
                     hi_params["centsub_sigma_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_sigma"],
               yerr=[co_params["centsub_sigma_low_lim"],
                     co_params["centsub_sigma_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(1.3, 14, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].errorbar(bin_cents, hi_params["peaksub_sigma"],
               yerr=[hi_params["peaksub_sigma_low_lim"],
                     hi_params["peaksub_sigma_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_sigma"],
               yerr=[co_params["peaksub_sigma_low_lim"],
                     co_params["peaksub_sigma_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
# ax[2].set_xlabel("Radius (kpc)")
ax[2].text(1.3, 14, "Peak Vel. subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
p.tight_layout()
p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.savefig(osjoin(figure_folder, "total_profile_radial_widths_HI_CO21.pdf"))
fig.savefig(osjoin(figure_folder, "total_profile_radial_widths_HI_CO21.png"))

p.close()

# Peak velocity only
onecolumn_figure()

p.errorbar(bin_cents, hi_params["peaksub_sigma"],
           yerr=[hi_params["peaksub_sigma_low_lim"],
                 hi_params["peaksub_sigma_up_lim"]],
           label='HI', color=cpal[2],
           drawstyle='steps-mid')
p.errorbar(bin_cents, co_params["peaksub_sigma"],
           yerr=[co_params["peaksub_sigma_low_lim"],
                 co_params["peaksub_sigma_up_lim"]],
           linestyle='--', label='CO(2-1)',
           drawstyle='steps-mid', color=cpal[1])
p.legend(loc='lower right', frameon=True)
p.grid()
p.ylim([0, 9])
p.xlabel("Radius (kpc)")
p.ylabel(r"$\sigma_{\rm HWHM}$ (km/s)")

p.tight_layout()

p.savefig(osjoin(figure_folder, "total_profile_radial_widths_peakvel_HI_CO21.pdf"))
p.savefig(osjoin(figure_folder, "total_profile_radial_widths_peakvel_HI_CO21.png"))

p.close()

onecolumn_Npanel_figure(N=2)

fig, axs = p.subplots(2, 1, sharex=True)

axs[0].errorbar(bin_cents, hi_params["peaksub_sigma"],
                yerr=[hi_params["peaksub_sigma_low_lim"],
                      hi_params["peaksub_sigma_up_lim"]],
                label='HI', color=cpal[2],
                drawstyle='steps-mid')
axs[0].errorbar(bin_cents, co_params["peaksub_sigma"],
                yerr=[co_params["peaksub_sigma_low_lim"],
                      co_params["peaksub_sigma_up_lim"]],
                linestyle='--', label='CO(2-1)',
                drawstyle='steps-mid', color=cpal[1])
axs[0].legend(loc='lower right', frameon=True)
axs[0].grid()
axs[0].set_ylim([0, 9])
axs[0].set_ylabel(r"$\sigma_{\rm HWHM}$ (km/s)")

axs[1].errorbar(bin_cents, hi_params["peaksub_f_wings"],
                yerr=[hi_params["peaksub_f_wings_low_lim"],
                      hi_params["peaksub_f_wings_up_lim"]],
                label='HI', color=cpal[2],
                drawstyle='steps-mid')
axs[1].errorbar(bin_cents, co_params["peaksub_f_wings"],
                yerr=[co_params["peaksub_f_wings_low_lim"],
                      co_params["peaksub_f_wings_up_lim"]],
                linestyle='--', label='CO(2-1)',
                drawstyle='steps-mid', color=cpal[1])
axs[1].grid()
axs[1].set_ylim([-0.1, 0.5])
axs[1].set_ylabel(r"$f_{\rm wings}$")

axs[1].set_xlabel("Radius (kpc)")


p.tight_layout()

p.savefig(osjoin(figure_folder, "total_profile_radial_widths_fwings_peakvel_HI_CO21.pdf"))
p.savefig(osjoin(figure_folder, "total_profile_radial_widths_fwings_peakvel_HI_CO21.png"))

p.close()


# Compare velocity at the peak of the stacks
twocolumn_twopanel_figure()

fig, ax = p.subplots(1, 3, sharey=True)

ax[0].errorbar(bin_cents, hi_params["rotsub_v_peak"],
               yerr=[hi_params["rotsub_v_peak_low_lim"],
                     hi_params["rotsub_v_peak_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_v_peak"],
               yerr=[co_params["rotsub_v_peak_low_lim"],
                     co_params["rotsub_v_peak_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].grid()
# ax[0].set_ylim([0.25, 16])
ax[0].text(1.3, 10, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel("Central Velocity (km/s)")

ax[1].errorbar(bin_cents, hi_params["centsub_v_peak"],
               yerr=[hi_params["centsub_v_peak_low_lim"],
                     hi_params["centsub_v_peak_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_v_peak"],
               yerr=[co_params["centsub_v_peak_low_lim"],
                     co_params["centsub_v_peak_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(1.3, 10, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].errorbar(bin_cents, hi_params["peaksub_v_peak"],
               yerr=[hi_params["peaksub_v_peak_low_lim"],
                     hi_params["peaksub_v_peak_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_v_peak"],
               yerr=[co_params["peaksub_v_peak_low_lim"],
                     co_params["peaksub_v_peak_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
ax[2].legend(loc='lower left', frameon=True)
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
               yerr=[hi_params["rotsub_f_wings_low_lim"],
                     hi_params["rotsub_f_wings_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_f_wings"],
               yerr=[co_params["rotsub_f_wings_low_lim"],
                     co_params["rotsub_f_wings_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].grid()
ax[0].set_ylim([-0.33, 0.75])
ax[0].text(1.3, 0.58, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
# ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel(r"$f_{\rm wings}$")
ax[0].axhline(0., color='k', linestyle=':')

ax[1].errorbar(bin_cents, hi_params["centsub_f_wings"],
               yerr=[hi_params["centsub_f_wings_low_lim"],
                     hi_params["centsub_f_wings_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_f_wings"],
               yerr=[co_params["centsub_f_wings_low_lim"],
                     co_params["centsub_f_wings_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(1.3, 0.58, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[1].axhline(0., color='k', linestyle=':')

ax[2].errorbar(bin_cents, hi_params["peaksub_f_wings"],
               yerr=[hi_params["peaksub_f_wings_low_lim"],
                     hi_params["peaksub_f_wings_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_f_wings"],
               yerr=[co_params["peaksub_f_wings_low_lim"],
                     co_params["peaksub_f_wings_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
ax[2].legend(loc='lower left', frameon=True)
# ax[2].set_xlabel("Radius (kpc)")
ax[2].text(1.3, 0.58, "Peak Vel. subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[2].axhline(0., color='k', linestyle=':')
p.tight_layout()
p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.savefig(osjoin(figure_folder, "total_profile_radial_fwings_HI_CO21.pdf"))
fig.savefig(osjoin(figure_folder, "total_profile_radial_fwings_HI_CO21.png"))

p.close()

# Asymmetry: this is funky. The profiles are genuinely asymmetrical, and the
# Druard stacked profiles show the same thing (I checked by-eye)
# What does that say? Is the \co peak velocity farther from the \hi at larger
# radii? But why would it only be in one direction?
# The asymmetry is consistent, but not significant compared to the uncertainty

twocolumn_figure(font_scale=1.3)
fig, ax = p.subplots(1, 3, sharey=True)

ax[0].errorbar(bin_cents, hi_params["rotsub_asymm"],
               yerr=[hi_params["rotsub_asymm_low_lim"],
                     hi_params["rotsub_asymm_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_asymm"],
               yerr=[co_params["rotsub_asymm_low_lim"],
                     co_params["rotsub_asymm_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].grid()
ax[0].set_ylim([-0.7, 1.1])
ax[0].text(0.5, 0.8, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
# ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel(r"Asymmetry")

ax[1].errorbar(bin_cents, hi_params["centsub_asymm"],
               yerr=[hi_params["centsub_asymm_low_lim"],
                     hi_params["centsub_asymm_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_asymm"],
               yerr=[co_params["centsub_asymm_low_lim"],
                     co_params["centsub_asymm_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(0.5, 0.8, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].errorbar(bin_cents, hi_params["peaksub_asymm"],
               yerr=[hi_params["peaksub_asymm_low_lim"],
                     hi_params["peaksub_asymm_up_lim"]],
               label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_asymm"],
               yerr=[co_params["peaksub_asymm_low_lim"],
                     co_params["peaksub_asymm_up_lim"]],
               linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
ax[2].legend(loc='lower left', frameon=True)
# ax[2].set_xlabel("Radius (kpc)")
ax[2].text(0.5, 0.8, "Peak Vel. subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
p.tight_layout()
p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.savefig(osjoin(figure_folder, "total_profile_radial_asymm_HI_CO21.pdf"))
fig.savefig(osjoin(figure_folder, "total_profile_radial_asymm_HI_CO21.png"))

p.close()

# The other parameters don't have interesting features.
# sigma_wing is highly uncertain.
# kappa are all consistent with a Gaussian shape.

# Compare the N/S profiles to each other.

twocolumn_figure(font_scale=1.3)

ylims = [[2, 15], None, [-1, 1], None, [-0.5, 0.5], [-0.2, 0.2]]

for ctr, param in enumerate(param_names):
    p.figure(1, figsize=(8.4, 11)).clf()

    fig, ax = p.subplots(3, 3,
                         sharex=True,
                         sharey=True, num=1)

    ax[0, 0].errorbar(bin_cents, hi_params["rotsub_{}".format(param)],
                      yerr=[hi_params["rotsub_{}_low_lim".format(param)],
                            hi_params["rotsub_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[0, 0].errorbar(bin_cents, co_params["rotsub_{}".format(param)],
                      yerr=[co_params["rotsub_{}_low_lim".format(param)],
                            co_params["rotsub_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[0, 0].grid()
    if ylims[ctr] is not None:
        ax[0, 0].set_ylim(ylims[ctr])
    # ax[0, 0].set_xlabel("Radius (kpc)")
    ax[0, 0].set_ylabel(r"Rot. Sub. {}".format(param.replace("_", " ")))

    ax[0, 0].set_xlabel('Whole')
    ax[0, 0].xaxis.set_label_position('top')

    ax[0, 1].errorbar(bin_cents, hi_params["rotsub_n_{}".format(param)],
                      yerr=[hi_params["rotsub_n_{}_low_lim".format(param)],
                            hi_params["rotsub_n_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[0, 1].errorbar(bin_cents, co_params["rotsub_n_{}".format(param)],
                      yerr=[co_params["rotsub_n_{}_low_lim".format(param)],
                            co_params["rotsub_n_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[0, 1].grid()
    ax[0, 1].set_xlabel('North')
    ax[0, 1].xaxis.set_label_position('top')

    ax[0, 2].errorbar(bin_cents, hi_params["rotsub_s_{}".format(param)],
                      yerr=[hi_params["rotsub_s_{}_low_lim".format(param)],
                            hi_params["rotsub_s_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[0, 2].errorbar(bin_cents, co_params["rotsub_s_{}".format(param)],
                      yerr=[co_params["rotsub_s_{}_low_lim".format(param)],
                            co_params["rotsub_s_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[0, 2].grid()
    ax[0, 2].set_xlabel('South')
    ax[0, 2].xaxis.set_label_position('top')

    ax[1, 0].errorbar(bin_cents, hi_params["centsub_{}".format(param)],
                      yerr=[hi_params["centsub_{}_low_lim".format(param)],
                            hi_params["centsub_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[1, 0].errorbar(bin_cents, co_params["centsub_{}".format(param)],
                      yerr=[co_params["centsub_{}_low_lim".format(param)],
                            co_params["centsub_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[1, 0].grid()
    ax[1, 0].set_ylabel(r"Cent. Sub. {}".format(param.replace("_", " ")))

    ax[1, 1].errorbar(bin_cents, hi_params["centsub_n_{}".format(param)],
                      yerr=[hi_params["centsub_n_{}_low_lim".format(param)],
                            hi_params["centsub_n_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[1, 1].errorbar(bin_cents, co_params["centsub_n_{}".format(param)],
                      yerr=[co_params["centsub_n_{}_low_lim".format(param)],
                            co_params["centsub_n_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[1, 1].grid()

    ax[1, 2].errorbar(bin_cents, hi_params["centsub_s_{}".format(param)],
                      yerr=[hi_params["centsub_s_{}_low_lim".format(param)],
                            hi_params["centsub_s_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[1, 2].errorbar(bin_cents, co_params["centsub_s_{}".format(param)],
                      yerr=[co_params["centsub_s_{}_low_lim".format(param)],
                            co_params["centsub_s_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[1, 2].grid()

    ax[2, 0].errorbar(bin_cents, hi_params["peaksub_{}".format(param)],
                      yerr=[hi_params["peaksub_{}_low_lim".format(param)],
                            hi_params["peaksub_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[2, 0].errorbar(bin_cents, co_params["peaksub_{}".format(param)],
                      yerr=[co_params["peaksub_{}_low_lim".format(param)],
                            co_params["peaksub_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[2, 0].grid()
    ax[2, 0].legend(loc='upper left', frameon=True)
    # ax[2, 0].set_xlabel("Radius (kpc)")
    ax[2, 0].set_ylabel(r"Peak Sub. {}".format(param.replace("_", " ")))

    ax[2, 1].errorbar(bin_cents, hi_params["peaksub_n_{}".format(param)],
                      yerr=[hi_params["peaksub_n_{}_low_lim".format(param)],
                            hi_params["peaksub_n_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[2, 1].errorbar(bin_cents, co_params["peaksub_n_{}".format(param)],
                      yerr=[co_params["peaksub_n_{}_low_lim".format(param)],
                            co_params["peaksub_n_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[2, 1].grid()
    ax[2, 1].set_xlabel("Radius (kpc)")

    ax[2, 2].errorbar(bin_cents, hi_params["peaksub_s_{}".format(param)],
                      yerr=[hi_params["peaksub_s_{}_low_lim".format(param)],
                            hi_params["peaksub_s_{}_up_lim".format(param)]],
                      label='HI',
                      drawstyle='steps-mid')
    ax[2, 2].errorbar(bin_cents, co_params["peaksub_s_{}".format(param)],
                      yerr=[co_params["peaksub_s_{}_low_lim".format(param)],
                            co_params["peaksub_s_{}_up_lim".format(param)]],
                      linestyle='--', label='CO(2-1)',
                      drawstyle='steps-mid')
    ax[2, 2].grid()

    p.tight_layout()
    p.subplots_adjust(hspace=0.05,
                      wspace=0.05)

    # p.draw()
    # raw_input(param)

    if param == "sigma_wing":
        p.close()
        continue

    fig.savefig(osjoin(figure_folder, "total_profile_radial_N_S_{}_HI_CO21.pdf".format(param)))
    fig.savefig(osjoin(figure_folder, "total_profile_radial_N_S_{}_HI_CO21.png".format(param)))

    p.close()

# Compare the amount of excess line wing flux to the error beam estimate from
# Druard+14

error_beam_flux = (2.5e7 * u.solMass) / ((24 * u.pc).to(u.cm))**2 / \
    const.m_p.to(u.solMass) / (5e20 * u.cm**-2 / (u. K * u.km / u.s))

spec = OneDSpectrum.from_hdu(fits.open(co_stackpath("peakvel_stacked_{}.fits"
                                                    .format(maxrad_string))))

vel_diff = np.abs(np.diff(spec.spectral_axis[:2])[0]).to(u.km / u.s)

hwhm_mod = fit_hwhm(spec.spectral_axis.value, spec.value)

# Only include the wings within the velocity mask that was used for fitting
line_wing_flux = ((spec - hwhm_mod[-1](spec.spectral_axis.value) * u.K)[vel_mask].sum() *
                  vel_diff)

print("Error beam flux / Line wing excess flux: {}"
      .format(error_beam_flux / line_wing_flux))
# Error beam flux / Line wing excess flux: 0.449239614814

default_figure()
