
import astropy.units as u
from spectral_cube import SpectralCube, OneDSpectrum
from spectral_cube.cube_utils import largest_beam
import numpy as np
import matplotlib.pyplot as p
from pandas import DataFrame
from astropy.modeling import models, fitting
from astropy.io import fits
from os.path import join as osjoin
import os
import seaborn as sb


from cube_analysis.spectral_stacking import total_profile
from cube_analysis.spectral_stacking_models import fit_2gaussian, fit_hwhm

from paths import (fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, alltables_path, fourteenB_HI_data_path,
                   fourteenB_HI_data_wGBT_path)
from constants import hi_freq
from galaxy_params import gal_feath as gal
from plotting_styles import *

'''
Create profiles of HI after subtracting velocity surfaces.
'''

noSD_stackpath = lambda x: osjoin(fourteenB_HI_data_path("", no_check=True), "stacked_spectra", x)
comb_stackpath = lambda x: osjoin(fourteenB_HI_data_wGBT_path("", no_check=True), "stacked_spectra", x)


figure_folder = "stacked_profiles"
if not os.path.exists(allfigs_path(figure_folder)):
    os.mkdir(allfigs_path(figure_folder))

dr = 100 * u.pc
max_radius = (8.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)

nbins = np.int(np.floor(max_radius / dr))
inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)

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
    np.array([29., 86., 142., 200., 265., 321., 374., 430.,
              494., 543., 603., 670., 726., 782., 836., 890.,
              953., 1007., 1071., 1131., 1188., 1240., 1295., 1353.,
              1420., 1478., 1531., 1585., 1650., 1703., 1753., 1815.,
              1885., 1940., 1991., 2054., 2110., 2160., 2221., 2268.,
              2354., 2409., 2448., 2509., 2573., 2622., 2679., 2755.,
              2809., 2859., 2903., 2968., 3039., 3083., 3158., 3212.,
              3270., 3309., 3376., 3435., 3496., 3554., 3620., 3664.,
              3729., 3781., 3831., 3890., 3967., 4024., 4072., 4122.,
              4192., 4237., 4290., 4366., 4434., 4485., 4527., 4585.])
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

# Avg. noise level in the cube
sigma_noise = 2.8  # K


cpal = sb.color_palette()

# Load in the total profiles
total_spectrum_hi = OneDSpectrum.from_hdu(fits.open(noSD_stackpath("rotation_stacked.fits")))
total_spectrum_hi_cent = OneDSpectrum.from_hdu(fits.open(noSD_stackpath("centroid_stacked.fits")))
total_spectrum_hi_peakvel = OneDSpectrum.from_hdu(fits.open(noSD_stackpath("peakvel_stacked.fits")))

total_spectrum_hi_wGBT = OneDSpectrum.from_hdu(fits.open(comb_stackpath("rotation_stacked.fits")))
total_spectrum_hi_cent_wGBT = OneDSpectrum.from_hdu(fits.open(comb_stackpath("centroid_stacked.fits")))
total_spectrum_hi_peakvel_wGBT = OneDSpectrum.from_hdu(fits.open(comb_stackpath("peakvel_stacked.fits")))


spectra = [total_spectrum_hi, total_spectrum_hi_cent,
           total_spectrum_hi_peakvel,
           total_spectrum_hi_wGBT, total_spectrum_hi_cent_wGBT,
           total_spectrum_hi_peakvel_wGBT]
labels = ["Rot. Sub.", "Cent. Sub.", "Peak Sub.",
          "Feath. Rot. Sub.", "Feath. Cent. Sub.", "Feath. Peak Sub."]
file_labels = ["rotsub", "centsub", "peaksub",
               "rotsub_feath", "centsub_feath", "peaksub_feath"]
hi_fit_vals = {}
hi_fit_hwhm_vals = {}
hwhm_models = dict.fromkeys(file_labels)

# Model the spectra

onecolumn_figure()

for spectrum, label, file_label in zip(spectra, labels, file_labels):

    norm_intens = (spectrum / np.nanmax(spectrum)).value

    vels = spectrum.spectral_axis.value / 1000.
    parvals, parerrs, cov, parnames, g_HI = fit_2gaussian(vels, norm_intens)

    hi_fit_vals[label + " Params"] = parvals
    hi_fit_vals[label + " Errors"] = parerrs

    # HWHM fitting
    parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_HI_hwhm = \
        fit_hwhm(vels, norm_intens,
                 sigma_noise=sigma_noise / np.nanmax(spectrum.value),
                 nbeams=(num_pix_feath_total / npix_beam if "feath" in
                         file_label else num_pix_total / npix_beam))

    hi_fit_hwhm_vals[label + " Params"] = parvals_hwhm
    hi_fit_hwhm_vals[label + " Errors"] = parerrs_hwhm
    hwhm_models[file_label] = g_HI_hwhm

    # Note that the statistical errors on the mean are too small.
    ax1 = p.subplot2grid((5, 1), (0, 0), rowspan=4)
    p.plot(vels, norm_intens, '-', drawstyle='steps-mid')
    p.plot(vels, g_HI(vels), ':', label="Total Fit")
    p.plot(vels, g_HI["None_0"](vels), '--', label="Narrow")
    p.plot(vels, g_HI["None_1"](vels), '-.', label="Wide")
    p.ylabel("HI Normalized Intensity")
    p.xticks([])
    p.xlim([-50, 50])
    p.legend(frameon=True)
    p.ylim([-0.1, 1.1])
    p.grid()

    ax2 = p.subplot2grid((5, 1), (4, 0))
    p.plot(vels, norm_intens - g_HI(vels), '-', drawstyle='steps-mid')
    p.hlines(0.0, -100, 100, color='k')
    p.grid()
    # p.ylim([-0.015, 0.015])
    p.xlim([-50, 50])
    p.xlabel("Velocity (km/s)")
    p.ylabel("Residuals")

    p.tight_layout()
    p.subplots_adjust(hspace=0)

    filename = "total_profile_corrected_velocity_{}_hi_2gauss_fit".format(file_label)
    p.savefig(allfigs_path(osjoin(figure_folder, filename + ".pdf")))
    p.savefig(allfigs_path(osjoin(figure_folder, filename + ".png")))
    p.close()
    # p.draw()
    # raw_input("Next plot?")
    # p.clf()

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

# Plot the total profiles w/ and w/o GBT added.

twocolumn_figure()
fig, ax = p.subplots(2, 3, sharey=True, sharex=True)

ax[0, 0].plot(total_spectrum_hi.spectral_axis.to(u.km / u.s).value,
              (total_spectrum_hi / total_spectrum_hi.max()).value,
              '-', drawstyle='steps-mid', label="VLA")
vels = total_spectrum_hi.spectral_axis.to(u.km / u.s).value
ax[0, 0].fill_between(vels, hwhm_models['rotsub'](vels),
                      (total_spectrum_hi / total_spectrum_hi.max()).value,
                      facecolor='gray',
                      alpha=0.5)
ax[0, 0].axvline(hwhm_models['rotsub'].mean + hwhm_models['rotsub'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[0, 0].axvline(hwhm_models['rotsub'].mean - hwhm_models['rotsub'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[0, 0].text(-45, 0.88, "Rotation\nsubtracted",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 0].text(29, 0.88, "VLA",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 0].set_ylim([-0.02, 1.1])
ax[0, 0].set_xlim([-50, 50])
ax[0, 0].grid()

ax[0, 1].plot(total_spectrum_hi_cent.spectral_axis.to(u.km / u.s).value,
              (total_spectrum_hi_cent / total_spectrum_hi_cent.max()).value,
              '-', drawstyle='steps-mid')
vels = total_spectrum_hi_cent.spectral_axis.to(u.km / u.s).value
ax[0, 1].fill_between(vels, hwhm_models['centsub'](vels),
                      (total_spectrum_hi_cent / total_spectrum_hi_cent.max()).value,
                      facecolor='gray',
                      alpha=0.5)
ax[0, 1].axvline(hwhm_models['centsub'].mean + hwhm_models['centsub'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[0, 1].axvline(hwhm_models['centsub'].mean - hwhm_models['centsub'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[0, 1].text(-45, 0.88, "Centroid\nsubtracted",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 1].text(25, 0.88, "VLA",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 1].set_ylim([-0.02, 1.1])
ax[0, 1].set_xlim([-50, 50])
ax[0, 1].grid()

ax[0, 2].plot(total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value,
              (total_spectrum_hi_peakvel / total_spectrum_hi_peakvel.max()).value,
              '-', drawstyle='steps-mid')
vels = total_spectrum_hi_peakvel.spectral_axis.to(u.km / u.s).value
ax[0, 2].fill_between(vels, hwhm_models['peaksub'](vels),
                      (total_spectrum_hi_peakvel / total_spectrum_hi_peakvel.max()).value,
                      facecolor='gray',
                      alpha=0.5)
ax[0, 2].axvline(hwhm_models['peaksub'].mean + hwhm_models['peaksub'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[0, 2].axvline(hwhm_models['peaksub'].mean - hwhm_models['peaksub'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[0, 2].text(-45, 0.88, "Peak Vel.\nsubtracted",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 2].text(25, 0.88, "VLA",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 2].set_ylim([-0.02, 1.1])
ax[0, 2].set_xlim([-50, 50])
ax[0, 2].legend(frameon=True)
ax[0, 2].grid()

ax[1, 0].plot(total_spectrum_hi_wGBT.spectral_axis.to(u.km / u.s).value,
              (total_spectrum_hi_wGBT / total_spectrum_hi_wGBT.max()).value,
              '-', drawstyle='steps-mid', label="VLA+GBT")
vels = total_spectrum_hi_wGBT.spectral_axis.to(u.km / u.s).value
ax[1, 0].fill_between(vels, hwhm_models['rotsub_feath'](vels),
                      (total_spectrum_hi_wGBT / total_spectrum_hi_wGBT.max()).value,
                      facecolor='gray',
                      alpha=0.5)
ax[1, 0].axvline(hwhm_models['rotsub_feath'].mean + hwhm_models['rotsub_feath'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[1, 0].axvline(hwhm_models['rotsub_feath'].mean - hwhm_models['rotsub_feath'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[1, 0].text(-45, 0.88, "Rotation\nsubtracted",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 0].text(15, 0.88, "VLA+GBT",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 0].set_ylim([-0.02, 1.1])
ax[1, 0].set_xlim([-50, 50])
ax[1, 0].grid()

ax[1, 1].plot(total_spectrum_hi_cent_wGBT.spectral_axis.to(u.km / u.s).value,
              (total_spectrum_hi_cent_wGBT / total_spectrum_hi_cent_wGBT.max()).value,
              '-', drawstyle='steps-mid')
vels = total_spectrum_hi_cent_wGBT.spectral_axis.to(u.km / u.s).value
ax[1, 1].fill_between(vels, hwhm_models['centsub_feath'](vels),
                      (total_spectrum_hi_cent_wGBT / total_spectrum_hi_cent_wGBT.max()).value,
                      facecolor='gray',
                      alpha=0.5)
ax[1, 1].axvline(hwhm_models['centsub_feath'].mean + hwhm_models['centsub_feath'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[1, 1].axvline(hwhm_models['centsub_feath'].mean - hwhm_models['centsub_feath'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[1, 1].text(-45, 0.88, "Centroid\nsubtracted",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 1].text(15, 0.88, "VLA+GBT",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 1].set_xlabel("Velocity (km/s)")
ax[1, 1].set_ylim([-0.02, 1.1])
ax[1, 1].set_xlim([-50, 50])
ax[1, 1].grid()

ax[1, 2].plot(total_spectrum_hi_peakvel_wGBT.spectral_axis.to(u.km / u.s).value,
              (total_spectrum_hi_peakvel_wGBT / total_spectrum_hi_peakvel_wGBT.max()).value,
              '-', drawstyle='steps-mid')
vels = total_spectrum_hi_peakvel_wGBT.spectral_axis.to(u.km / u.s).value
ax[1, 2].fill_between(vels, hwhm_models['peaksub_feath'](vels),
                      (total_spectrum_hi_peakvel_wGBT / total_spectrum_hi_peakvel_wGBT.max()).value,
                      facecolor='gray',
                      alpha=0.5)
ax[1, 2].axvline(hwhm_models['peaksub_feath'].mean + hwhm_models['peaksub_feath'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[1, 2].axvline(hwhm_models['peaksub_feath'].mean - hwhm_models['peaksub_feath'].stddev *
                 np.sqrt(2 * np.log(2)), linestyle='--', color=cpal[1])
ax[1, 2].text(-45, 0.88, "Peak Vel.\nsubtracted",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 2].text(15, 0.88, "VLA+GBT",
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 2].set_ylim([-0.02, 1.1])
ax[1, 2].set_xlim([-50, 50])
ax[1, 2].legend(frameon=True)
ax[1, 2].grid()

# Add a ylabel
fig.text(0.01, 0.5, 'Normalized Total Intensity',
         va='center', rotation='vertical')

p.tight_layout()
p.subplots_adjust(hspace=0.02,
                  wspace=0.02,
                  left=0.08)

fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_corrected_velocity_HI.pdf")))
fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_corrected_velocity_HI.png")))
p.close()

# Save parameter table
hi_param_df = DataFrame(hi_fit_vals, index=parnames)
hi_param_df.to_latex(alltables_path("hi_gaussian_totalprof_fits.tex"))
hi_param_df.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits.csv",
                                          no_check=True))

hi_hwhm_param_df = DataFrame(hi_fit_hwhm_vals, index=parnames_hwhm)
hi_hwhm_param_df.to_latex(alltables_path("hi_gaussian_totalprof_hwhm_fits.tex"))
hi_hwhm_param_df.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_hwhm_fits.csv",
                                               no_check=True))


rot_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_north_{}.fits".format(wstring)))
rot_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_south_{}.fits".format(wstring)))
rot_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_{}.fits".format(wstring)))

cent_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_north_{}.fits".format(wstring)))
cent_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_south_{}.fits".format(wstring)))
cent_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_{}.fits".format(wstring)))

peakvel_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring)))
peakvel_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring)))
peakvel_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring)))

rot_stack_feath_n = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_north_{}.fits".format(wstring)))
rot_stack_feath_s = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_south_{}.fits".format(wstring)))
rot_stack_feath = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_{}.fits".format(wstring)))

cent_stack_feath_n = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_north_{}.fits".format(wstring)))
cent_stack_feath_s = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_south_{}.fits".format(wstring)))
cent_stack_feath = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_{}.fits".format(wstring)))

peakvel_stack_feath_n = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring)))
peakvel_stack_feath_s = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring)))
peakvel_stack_feath = \
    SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring)))

# How do the model parameters change with radius?

file_labels = ["rotsub", "rotsub_n", "rotsub_s",
               "centsub", "centsub_n", "centsub_s",
               "peaksub", "peaksub_n", "peaksub_s",
               "rotsub_feath", "rotsub_feath_n", "rotsub_feath_s",
               "centsub_feath", "centsub_feath_n", "centsub_feath_s",
               "peaksub_feath", "peaksub_feath_n", "peaksub_feath_s"]

hi_params = {}

param_names = parnames_hwhm

for sub in file_labels:
    for name in param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_error = "{}_stderr".format(par_name)

        hi_params[par_name] = np.zeros_like(inneredge.value)
        hi_params[par_error] = np.zeros_like(inneredge.value)


for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    hi_spectra = [rot_stack[:, ctr, 0],
                  rot_stack_n[:, ctr, 0],
                  rot_stack_s[:, ctr, 0],
                  cent_stack[:, ctr, 0],
                  cent_stack_n[:, ctr, 0],
                  cent_stack_s[:, ctr, 0],
                  peakvel_stack[:, ctr, 0],
                  peakvel_stack_n[:, ctr, 0],
                  peakvel_stack_s[:, ctr, 0],
                  rot_stack_feath[:, ctr, 0],
                  rot_stack_feath_n[:, ctr, 0],
                  rot_stack_feath_s[:, ctr, 0],
                  cent_stack_feath[:, ctr, 0],
                  cent_stack_feath_n[:, ctr, 0],
                  cent_stack_feath_s[:, ctr, 0],
                  peakvel_stack_feath[:, ctr, 0],
                  peakvel_stack_feath_n[:, ctr, 0],
                  peakvel_stack_feath_s[:, ctr, 0]]

    for spectrum, label in zip(hi_spectra, file_labels):

        vels = spectrum.spectral_axis.to(u.km / u.s).value

        if "_n" in label:
            nbeams = n_numpix[ctr] / npix_beam
        elif "_s" in label:
            nbeams = s_numpix[ctr] / npix_beam
        else:
            nbeams = num_pix[ctr] / npix_beam

        parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_HI_hwhm = \
            fit_hwhm(vels, spectrum.value, sigma_noise=sigma_noise,
                     nbeams=nbeams)

        for idx, name in enumerate(parnames_hwhm):
            par_name = "{0}_{1}".format(label, name)
            hi_params[par_name][ctr] = parvals_hwhm[idx]
            hi_params["{}_stderr".format(par_name)][ctr] = \
                parerrs_hwhm[idx]
            #     np.sqrt(cov[idx - idx_corr, idx - idx_corr])

bin_names = ["{}-{}".format(r0.value, r1)
             for r0, r1 in zip(inneredge, outeredge)]

hi_radial_fits = DataFrame(hi_params, index=bin_names)

hi_radial_fits.to_latex(alltables_path("hi_gaussian_hwhm_totalprof_fits_radial.tex"))
hi_radial_fits.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_hwhm_totalprof_fits_radial.csv",
                                             no_check=True))

# Plot comparisons of these fits -- showing f_wings, kappa, asymm
#
bin_cents = (outeredge - dr / 2.).to(u.kpc).value

twocolumn_figure()

fig, ax = p.subplots(2, 3, sharex=True)

ax[0, 0].errorbar(bin_cents, hi_params["rotsub_f_wings"],
                  yerr=hi_params["rotsub_f_wings_stderr"],
                  color=cpal[0], label='Rot. Sub.',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[0, 0].errorbar(bin_cents, hi_params["centsub_f_wings"],
                  yerr=hi_params["centsub_f_wings_stderr"],
                  color=cpal[1], label='Cent. Sub.',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[0, 0].errorbar(bin_cents, hi_params["peaksub_f_wings"],
                  yerr=hi_params["peaksub_f_wings_stderr"],
                  color=cpal[2], label='Peak Sub.',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[0, 0].fill_between([0, 0.5], -0.1, 0.37, facecolor='gray', alpha=0.5)
ax[0, 0].set_ylim([-0.1, 0.37])
ax[0, 0].set_xlim([0, 8])
ax[0, 0].grid()
ax[0, 0].text(3.3, -0.06, "VLA",
              bbox={"boxstyle": "square", "facecolor": "w"})

# ax[0].set_ylim([0.25, 15])
ax[0, 0].set_ylabel(r"f$_{\rm wings}$")


ax[1, 0].errorbar(bin_cents, hi_params["rotsub_feath_f_wings"],
                  yerr=hi_params["rotsub_feath_f_wings_stderr"],
                  color=cpal[0], label='Rot. Sub.',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[1, 0].errorbar(bin_cents, hi_params["centsub_feath_f_wings"],
                  yerr=hi_params["centsub_feath_f_wings_stderr"],
                  color=cpal[1], label='Cent. Sub.',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[1, 0].errorbar(bin_cents, hi_params["peaksub_feath_f_wings"],
                  yerr=hi_params["peaksub_feath_f_wings_stderr"],
                  color=cpal[2], label='Peak Sub.',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[1, 0].grid()
ax[1, 0].set_ylabel(r"f$_{\rm wings}$")
ax[1, 0].fill_between([0, 0.5], -0.1, 0.37, facecolor='gray', alpha=0.5)
ax[1, 0].set_ylim([-0.1, 0.37])
ax[1, 0].text(2.2, -0.06, "VLA + GBT",
              bbox={"boxstyle": "square", "facecolor": "w"})

# kappa
ax[0, 1].errorbar(bin_cents, hi_params["rotsub_kappa"],
                  yerr=hi_params["rotsub_kappa_stderr"],
                  color=cpal[0], label='Rot. Sub.',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[0, 1].errorbar(bin_cents, hi_params["centsub_kappa"],
                  yerr=hi_params["centsub_kappa_stderr"],
                  color=cpal[1], label='Cent. Sub.',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[0, 1].errorbar(bin_cents, hi_params["peaksub_kappa"],
                  yerr=hi_params["peaksub_kappa_stderr"],
                  color=cpal[2], label='Peak Sub.',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[0, 1].grid()
ax[0, 1].fill_between([0, 0.5], -0.09, 0.03, facecolor='gray', alpha=0.5)
ax[0, 1].set_ylim([-0.09, 0.03])

# ax[0].set_ylim([0.25, 15])
ax[0, 1].set_ylabel(r"$\kappa$")


ax[1, 1].errorbar(bin_cents, hi_params["rotsub_feath_kappa"],
                  yerr=hi_params["rotsub_feath_kappa_stderr"],
                  color=cpal[0], label='Rot. Sub.',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[1, 1].errorbar(bin_cents, hi_params["centsub_feath_kappa"],
                  yerr=hi_params["centsub_feath_kappa_stderr"],
                  color=cpal[1], label='Cent. Sub.',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[1, 1].errorbar(bin_cents, hi_params["peaksub_feath_kappa"],
                  yerr=hi_params["peaksub_feath_kappa_stderr"],
                  color=cpal[2], label='Peak Sub.',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[1, 1].grid()
ax[1, 1].set_ylabel(r"$\kappa$")
ax[1, 1].set_xlabel("Radius (kpc)")
ax[1, 1].fill_between([0, 0.5], -0.09, 0.03, facecolor='gray', alpha=0.5)
ax[1, 1].set_ylim([-0.09, 0.03])

# kappa
ax[0, 2].errorbar(bin_cents, hi_params["rotsub_asymm"],
                  yerr=hi_params["rotsub_asymm_stderr"],
                  color=cpal[0], label='Rot. Sub.',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[0, 2].errorbar(bin_cents, hi_params["centsub_asymm"],
                  yerr=hi_params["centsub_asymm_stderr"],
                  color=cpal[1], label='Cent. Sub.',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[0, 2].errorbar(bin_cents, hi_params["peaksub_asymm"],
                  yerr=hi_params["peaksub_asymm_stderr"],
                  color=cpal[2], label='Peak Sub.',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[0, 2].grid()
ax[0, 2].legend(frameon=True)
# ax[0].set_ylim([0.25, 15])
ax[0, 2].set_ylabel(r"Asymmetry")
ax[0, 2].fill_between([0, 0.5], -0.01, 0.36, facecolor='gray', alpha=0.5)
ax[0, 2].set_ylim([-0.01, 0.36])


ax[1, 2].errorbar(bin_cents, hi_params["rotsub_feath_asymm"],
                  yerr=hi_params["rotsub_feath_asymm_stderr"],
                  color=cpal[0], label='Rot. Sub.',
                  linestyle='-',
                  drawstyle='steps-mid')
ax[1, 2].errorbar(bin_cents, hi_params["centsub_feath_asymm"],
                  yerr=hi_params["centsub_feath_asymm_stderr"],
                  color=cpal[1], label='Cent. Sub.',
                  linestyle='--',
                  drawstyle='steps-mid')
ax[1, 2].errorbar(bin_cents, hi_params["peaksub_feath_asymm"],
                  yerr=hi_params["peaksub_feath_asymm_stderr"],
                  color=cpal[2], label='Peak Sub.',
                  linestyle='-.',
                  drawstyle='steps-mid')
ax[1, 2].grid()
ax[1, 2].set_ylabel(r"Asymmetry")
ax[1, 2].fill_between([0, 0.5], -0.01, 0.36, facecolor='gray', alpha=0.5)
ax[1, 2].set_ylim([-0.01, 0.36])

p.tight_layout()

fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_radial_params_HI.pdf")))
fig.savefig(allfigs_path(osjoin(figure_folder, "total_profile_radial_params_HI.png")))

p.close()

default_figure()
