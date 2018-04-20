import astropy.units as u
from spectral_cube import OneDSpectrum, SpectralCube
import numpy as np
from pandas import DataFrame
from astropy.io import fits
from astropy.table import Table, Column
from os.path import join as osjoin
import os
from radio_beam import Beam
from astropy import log

from cube_analysis.spectral_stacking_models import fit_hwhm

from paths import (allfigs_path, alltables_path, fourteenB_HI_data_wGBT_path)

'''
Model and compare the stacked profiles from total_stacked_profiles_lowres.py
'''

figure_folder = allfigs_path("stacked_profiles")
if not os.path.exists(figure_folder):
    os.mkdir(figure_folder)


dr = 200 * u.pc
max_radius = (10.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)
maxrad_string = "{0}{1}".format(int(max_radius.value), max_radius.unit)

# Load the total HI profiles in
hi_stackpath = lambda x: osjoin(fourteenB_HI_data_wGBT_path("smooth_2beam/stacked_spectra"), x)

total_spectrum_hi_cent = OneDSpectrum.from_hdu(fits.open(hi_stackpath("centroid_stacked_{}.fits".format(maxrad_string))))
total_spectrum_hi_peakvel = OneDSpectrum.from_hdu(fits.open(hi_stackpath("peakvel_stacked_{}.fits".format(maxrad_string))))

hi_spectra = [total_spectrum_hi_cent,
              total_spectrum_hi_peakvel]
hi_fit_vals = {}

labels = ["Cent. Sub.", "Peak Sub."]
file_labels = ["centsub", "peaksub"]

sigma_noise_hi = 2.8  # K

# Pixels in each radial bin
num_pix = np.load(hi_stackpath("radial_stacking_pixelsinbin_{}.npy".format(wstring)))

num_pix_total = num_pix.sum()

npix_beam = float((Beam(37.983341217 * u.arcsec).as_tophat_kernel((3 * u.arcsec).to(u.deg)).array > 0).sum())
# 177 pixels

log.info("Modelling total spectra. 2 * beam")
for hi_spectrum, label, file_label in zip(hi_spectra, labels, file_labels):

    # HWHM fitting
    hi_vels = hi_spectrum.spectral_axis.to(u.km / u.s).value
    vel_mask_hi = np.logical_and(hi_vels < 50, hi_vels > -50)

    parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_HI = \
        fit_hwhm(hi_vels, hi_spectrum.value,
                 sigma_noise=sigma_noise_hi,
                 nbeams=num_pix_total / npix_beam,
                 niters=100)

    hi_fit_vals[label + " Params"] = parvals_hwhm
    hi_fit_vals[label + " Lower Limit"] = parerrs_hwhm[0]
    hi_fit_vals[label + " Upper Limit"] = parerrs_hwhm[1]

# Make tables folders if needed
hi_tab_path = fourteenB_HI_data_wGBT_path("smooth_2beam/tables", no_check=True)
if not os.path.exists(hi_tab_path):
    os.mkdir(hi_tab_path)

hi_param_df = DataFrame(hi_fit_vals, index=parnames_hwhm)
hi_param_df.to_latex(alltables_path("hi_gaussian_totalprof_fits_hwhm_38arcsec.tex"))
hi_param_df.to_csv(fourteenB_HI_data_wGBT_path("smooth_2beam/tables/hi_gaussian_totalprof_fits_hwhm_38arcsec_{}.csv".format(maxrad_string),
                                               no_check=True))

# Same for the radial bins.

# Load in radial profiles

total_spectrum_hi_radial_cent = SpectralCube.read(hi_stackpath("centroid_stacked_radial_{}.fits".format(wstring)))

total_spectrum_hi_radial_peakvel = SpectralCube.read(hi_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)))

labels = ["centsub", "peaksub"]

# Make the bin edges
# Convert these into kpc
inneredge = np.arange(0, max_radius.value, dr.value) / 1000.
outeredge = (inneredge + dr.to(u.kpc).value)

# How do the model parameters change with radius?
hi_params = {}

hi_models = {}

param_names = ["sigma", "v_peak", "f_wings", "sigma_wing", "asymm", "kappa"]

for sub in labels:
    for name in param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_lowlim = "{}_low_lim".format(par_name)
        par_uplim = "{}_up_lim".format(par_name)

        hi_params[par_name] = np.zeros_like(inneredge)
        hi_params[par_lowlim] = np.zeros_like(inneredge)
        hi_params[par_uplim] = np.zeros_like(inneredge)

    hi_models[sub] = []

log.info("Modelling radial spectra. 2 * beam")
for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):
    print(ctr, len(inneredge))
    hi_spectra = [total_spectrum_hi_radial_cent[:, ctr, 0],
                  total_spectrum_hi_radial_peakvel[:, ctr, 0]]

    pix_in_bin = num_pix[ctr]

    for spectrum, label in zip(hi_spectra, labels):

        vels = spectrum.spectral_axis.to(u.km / u.s).value
        norm_intens = (spectrum / spectrum.max()).value

        vel_mask = np.logical_and(vels < 50, vels > -50)

        params, stderrs, names, hwhm_mod = \
            fit_hwhm(vels[vel_mask], spectrum.value[vel_mask],
                     sigma_noise=sigma_noise_hi,
                     nbeams=pix_in_bin / npix_beam,
                     niters=100)
        hi_models[label].append(hwhm_mod)

        for idx, name in enumerate(names):
            par_name = "{0}_{1}".format(label, name)
            hi_params[par_name][ctr] = params[idx]
            hi_params["{}_low_lim".format(par_name)][ctr] = \
                np.abs(stderrs[0, idx])
            hi_params["{}_up_lim".format(par_name)][ctr] = \
                np.abs(stderrs[1, idx])

bin_names = ["{}-{}".format(r0, r1)
             for r0, r1 in zip(inneredge, outeredge)]

hi_radial_fits = DataFrame(hi_params, index=bin_names)

hi_radial_fits.to_latex(alltables_path("hi_hwhm_totalprof_fits_38arcsec_radial_{}.tex".format(wstring)))
hi_radial_fits.to_csv(fourteenB_HI_data_wGBT_path("smooth_2beam/tables/hi_hwhm_totalprof_fits_38arcsec_radial_{}.csv".format(wstring),
                                                  no_check=True))
