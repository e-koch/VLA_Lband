
'''
The H2/HI ratio depends critically on the velocity extent used for both.
CO is straightforward, but HI is not. This script explores a couple methods of
using the CO velocity extent as a prior for the HI.

-

'''

from spectral_cube import SpectralCube
import numpy as np
import pylab as pl
from astropy.io import fits
from astropy.convolution import Box1DKernel
from astropy.table import Table
import astropy.units as u
from astropy.utils.console import ProgressBar
import os
from os.path import join as osjoin
import seaborn as sb
from scipy.special import erf

from cube_analysis.spectral_stacking_models import find_hwhm, fit_gaussian

from paths import (fourteenB_wGBT_HI_file_dict, iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path, allfigs_path)
from plotting_styles import (default_figure, onecolumn_figure,
                             onecolumn_twopanel_figure)
from galaxy_params import gal_feath as gal
from constants import (co21_mass_conversion, hi_mass_conversion, hi_freq,
                       beam_eff_30m_druard)


default_figure()
cpal = sb.color_palette()

# cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['PeakSub_Cube'])
cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])
co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.fits"))
co_rms = fits.open(iram_co21_14B088_data_path('m33.rms.14B-088_HI.fits'))[0]


co_mask = fits.open(iram_co21_14B088_data_path(
    "m33.co21_iram.14B-088_HI_source_mask.fits"))[0].data

# Require there be 3 pixels in the mask
co_mask_spatial = co_mask.sum(0) > 3
yposns, xposns = np.where(co_mask_spatial)

hi_specaxis = cube.spectral_axis
co_specaxis = co_cube.spectral_axis

hi_chanwidth = np.abs(np.diff(hi_specaxis[:2])[0].value)
co_chanwidth = np.abs(np.diff(co_specaxis[:2])[0].value)

# Extend the HI fitting mask by 3 channels
mask_pad = hi_chanwidth * 3

# Smooth the HI prior to fitting. Smoothing across 10 channels ~ 2 km/s
kern = Box1DKernel(10)

plot_spectra = False

results = {"amp_CO": np.zeros((yposns.shape)) * u.K,
           "mean_CO": np.zeros((yposns.shape)) * u.m / u.s,
           "sigma_CO": np.zeros((yposns.shape)) * u.m / u.s,
           "amp_stderr_CO": np.zeros((yposns.shape)) * u.K,
           "mean_stderr_CO": np.zeros((yposns.shape)) * u.m / u.s,
           "sigma_stderr_CO": np.zeros((yposns.shape)) * u.m / u.s,
           "amp_HI": np.zeros((yposns.shape)) * u.K,
           "mean_HI": np.zeros((yposns.shape)) * u.m / u.s,
           "sigma_HI": np.zeros((yposns.shape)) * u.m / u.s,
           "amp_stderr_HI": np.zeros((yposns.shape)) * u.K,
           "mean_stderr_HI": np.zeros((yposns.shape)) * u.m / u.s,
           "sigma_stderr_HI": np.zeros((yposns.shape)) * u.m / u.s,
           "amp_HI_nomask": np.zeros((yposns.shape)) * u.K,
           "mean_HI_nomask": np.zeros((yposns.shape)) * u.m / u.s,
           "sigma_HI_nomask": np.zeros((yposns.shape)) * u.m / u.s,
           "amp_stderr_HI_nomask": np.zeros((yposns.shape)) * u.K,
           "mean_stderr_HI_nomask": np.zeros((yposns.shape)) * u.m / u.s,
           "sigma_stderr_HI_nomask": np.zeros((yposns.shape)) * u.m / u.s,
           "coldens_HI_FWHM": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_HI_FWHM_stderr": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_CO_FWHM": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_CO_FWHM_stderr": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_HI_gauss": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_HI_gauss_stderr": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_HI_gauss_nomask": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_HI_gauss_nomask_stderr": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_CO_gauss": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_CO_gauss_stderr": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "multicomp_flag_HI": np.zeros(yposns.shape, dtype=bool),
           "multicomp_flag_CO": np.zeros(yposns.shape, dtype=bool)}

# Correct for the disk inclincation
inc = np.cos(gal.inclination)

# 30 m beam efficiency
beam_eff = beam_eff_30m_druard

# Roughly constant for all pixels
hi_err = 2.8 * u.K

# Pixels to save plots of as examples
# Order: simple both, multi-comp CO, both multi-comp,
# HI good fit w/ multi-comp, Bad HI multi-comp,
# HI w/ and w/o HWHM mask
example_pix = [(454, 924), (788, 576), (954, 438), (938, 681), (774, 680),
               (984, 523)]
example_params = []

# Make a new figures folder for these
figure_folder = allfigs_path("co_vs_hi/example_spectra/")
if not os.path.exists(allfigs_path(figure_folder)):
    os.mkdir(allfigs_path(figure_folder))

# rand_ord = np.random.choice(np.arange(len(yposns)), size=len(yposns),
#                             replace=False)

# for i, (y, x) in enumerate(ProgressBar(zip(yposns[rand_ord],
#                                            xposns[rand_ord]))):
for i, (y, x) in enumerate(ProgressBar(zip(yposns, xposns))):

    co_spectrum = co_cube[:, y, x]

    # The CO error can vary. Pull from the RMS map
    co_err = co_rms.data[y, x] * u.K

    # First beam is the largest. Change is far far smaller than 1 pixel
    hi_spectrum = cube[:, y, x].to(u.K, cube.beams[0].jtok_equiv(hi_freq))

    # Try fitting with the spectral response included
    # From Sun+18
    r = 0.26
    k = 0.47 * r - 0.23 * r**2 - 0.16 * r**3 + 0.43 * r**4

    co_params, co_stderrs, co_cov, co_parnames, co_model = \
        fit_gaussian(co_specaxis, co_spectrum.quantity,
                     sigma=co_err, use_discrete=True,
                     kernel=np.array([k, 1 - 2 * k, k]),
                     add_chan_width_err=False)

    if np.isnan(co_cov).any():
        results["multicomp_flag_CO"][i] = True

    results["amp_CO"][i] = co_params[0] * u.K
    results["mean_CO"][i] = co_params[1] * u.m / u.s
    results["sigma_CO"][i] = co_params[2] * u.m / u.s
    # Adjust the mean and sigma to account for finite channel widths
    # Take to be half of the channel width.
    results["amp_stderr_CO"][i] = co_stderrs[0] * u.K
    results["mean_stderr_CO"][i] = co_stderrs[1] * u.m / u.s
    results["sigma_stderr_CO"][i] = co_stderrs[2] * u.m / u.s

    # Make a mask centered around the CO peak with 5x the CO FWHM
    # This will hopefully help centre the HI fit to the CO if there is
    # a brighter velocity component.
    co_hwhm = co_params[-1] * np.sqrt(2 * np.log(2))
    co_mask = np.logical_and(hi_specaxis.value >= co_params[1] - 3 * co_hwhm,
                             hi_specaxis.value <= co_params[1] + 3 * co_hwhm)

    # Limit the HI fitting to the peak. Smooth the HI spectrum to get a good
    # HWHM estimate
    # This will fail for spectra with multiple components blended together.
    # Flag those that fail and make some rough estimates instead.
    try:
        hi_sigest, hi_hwhm = find_hwhm(hi_specaxis[co_mask],
                                       hi_spectrum.spectral_smooth(kern)[co_mask])[:2]
    except ValueError:
        # Just adopt some typical values if the estimate fails. We'll assign
        # those failures to blended profiles
        hi_sigest = 7000.
        hwhm_factor = np.sqrt(2 * np.log(2))
        peak_posn = np.nanargmax(hi_spectrum.spectral_smooth(kern)[co_mask])
        hi_hwhm = np.array([hi_specaxis[co_mask][peak_posn].value - hi_sigest * hwhm_factor,
                            hi_specaxis[co_mask][peak_posn].value + hi_sigest * hwhm_factor])
        results["multicomp_flag_HI"][i] = True

    hwhm_mask = np.logical_and(hi_specaxis.value >= hi_hwhm[0] - mask_pad,
                               hi_specaxis.value <= hi_hwhm[1] + mask_pad)

    # Give some initial guesses
    p0 = (np.nanmax(hi_spectrum[hwhm_mask]).value,
          hi_specaxis[hwhm_mask][np.nanargmax(hi_spectrum[hwhm_mask])].value,
          hi_sigest)

    hi_params, hi_stderrs, hi_cov, hi_parnames, hi_model = \
        fit_gaussian(hi_specaxis[hwhm_mask],
                     hi_spectrum.quantity[hwhm_mask], p0=p0,
                     sigma=hi_err,
                     add_chan_width_err=False)

    if np.isnan(hi_cov).any():
        results["multicomp_flag_HI"][i] = True

    # Catch bad HI fits where the fitted peak falls outside of the mask
    max_mask = hi_specaxis.value[np.where(hwhm_mask)[0][0]]
    min_mask = hi_specaxis.value[np.where(hwhm_mask)[0][-1]]
    if hi_params[1] < min_mask or hi_params[1] > max_mask:
        results['multicomp_flag_HI'][i] = True

    # By-eye, HI profiles with 1 major component have line widths of ~8 km/s
    # Let's assume those above 12 are due to bad fits.
    if hi_params[2] > 12000:
        results['multicomp_flag_HI'][i] = True

    # Use a similar criterion for the CO fits too. They have widths of ~3 to
    # 5 km /s. Set the cut-off at 8 km/s
    if co_params[2] > 8000:
        results['multicomp_flag_CO'][i] = True

    results["amp_HI"][i] = hi_params[0] * u.K
    results["mean_HI"][i] = hi_params[1] * u.m / u.s
    results["sigma_HI"][i] = hi_params[2] * u.m / u.s
    # Adjust the mean and sigma to account for finite channel widths
    # Take to be half of the channel width.
    results["amp_stderr_HI"][i] = hi_stderrs[0] * u.K
    results["mean_stderr_HI"][i] = hi_stderrs[1] * u.m / u.s
    results["sigma_stderr_HI"][i] = hi_stderrs[2] * u.m / u.s

    # How different are the fits if we fit a Gaussian to the whole HI line?
    # Good for comparing to many previous works

    hi_params_nomask, hi_stderrs_nomask, hi_cov_nomask, hi_parnames, \
        hi_model_nomask = \
        fit_gaussian(hi_specaxis,
                     hi_spectrum.quantity, p0=p0, sigma=hi_err)

    results["amp_HI_nomask"][i] = hi_params_nomask[0] * u.K
    results["mean_HI_nomask"][i] = hi_params_nomask[1] * u.m / u.s
    results["sigma_HI_nomask"][i] = hi_params_nomask[2] * u.m / u.s
    # Adjust the mean and sigma to account for finite channel widths
    # Take to be half of the channel width.
    results["amp_stderr_HI_nomask"][i] = hi_stderrs_nomask[0] * u.K
    results["mean_stderr_HI_nomask"][i] = hi_stderrs_nomask[1] * u.m / u.s
    results["sigma_stderr_HI_nomask"][i] = hi_stderrs_nomask[2] * u.m / u.s

    # Finally calculate column densities from:
    # 1) within the FWHM of each spectra,
    # 2) from the integral of the gaussian models

    # 1) These are divided by erf(sqrt(ln(2))) to account for the area outside
    # the FWHM

    fwhm_area_factor = erf(np.sqrt(np.log(2)))
    hi_hwhm_model = hi_params[-1] * np.sqrt(2 * np.log(2))
    hi_fwhm_mask = \
        np.logical_and(hi_specaxis.value >= hi_params[1] - hi_hwhm_model,
                       hi_specaxis.value <= hi_params[1] + hi_hwhm_model)
    results["coldens_HI_FWHM"][i] = \
        np.nansum(hi_spectrum[hi_fwhm_mask]) * \
        (hi_chanwidth / 1000. * u.km / u.s) * inc * hi_mass_conversion \
        / fwhm_area_factor

    results['coldens_HI_FWHM_stderr'][i] = sum(hi_fwhm_mask) * hi_err * \
        (hi_chanwidth / 1000. * u.km / u.s) * inc * hi_mass_conversion \
        / fwhm_area_factor

    co_fwhm_mask = \
        np.logical_and(co_specaxis.value >= co_params[1] - co_hwhm,
                       co_specaxis.value <= co_params[1] + co_hwhm)
    results["coldens_CO_FWHM"][i] = \
        np.nansum(co_spectrum[co_fwhm_mask]) * \
        (co_chanwidth / 1000. * u.km / u.s) * \
        inc * co21_mass_conversion / beam_eff / fwhm_area_factor

    results['coldens_CO_FWHM_stderr'][i] = sum(co_fwhm_mask) * co_err * \
        (co_chanwidth / 1000. * u.km / u.s) * inc * co21_mass_conversion \
        / beam_eff / fwhm_area_factor

    # 2)

    results["coldens_HI_gauss"][i] = \
        (hi_params[0] * u.K) * np.sqrt(2 * np.pi) * \
        (hi_params[-1] / 1000. * u.km / u.s) * inc * hi_mass_conversion

    results['coldens_HI_gauss_stderr'][i] = results["coldens_HI_gauss"][i] * \
        np.sqrt((results['amp_stderr_HI'][i] / results['amp_HI'][i])**2 +
                (results['sigma_stderr_HI'][i] / results['sigma_HI'][i])**2)

    results["coldens_HI_gauss_nomask"][i] = \
        (hi_params_nomask[0] * u.K) * np.sqrt(2 * np.pi) * \
        (hi_params_nomask[-1] / 1000. * u.km / u.s) * inc * hi_mass_conversion

    results['coldens_HI_gauss_nomask_stderr'][i] = \
        results["coldens_HI_gauss_nomask"][i] * \
        np.sqrt((results['amp_stderr_HI_nomask'][i] / results['amp_HI_nomask'][i])**2 +
                (results['sigma_stderr_HI_nomask'][i] / results['sigma_HI_nomask'][i])**2)

    results["coldens_CO_gauss"][i] = \
        (co_params[0] * u.K) * np.sqrt(2 * np.pi) * \
        (co_params[-1] / 1000. * u.km / u.s) * \
        inc * co21_mass_conversion / beam_eff

    results['coldens_CO_gauss_stderr'][i] = results["coldens_CO_gauss"][i] * \
        np.sqrt((results['amp_stderr_CO'][i] / results['amp_CO'][i])**2 +
                (results['sigma_stderr_CO'][i] / results['sigma_CO'][i])**2)

    if plot_spectra or (y, x) in example_pix:

        onecolumn_twopanel_figure()

        fig, ax = pl.subplots(3, 1, sharex=True, )

        vel_mask = np.logical_and(hi_specaxis.value >=
                                  max(np.min(hi_specaxis.value), hi_params[1] - 50000),
                                  hi_specaxis.value <=
                                  min(np.max(hi_specaxis.value), hi_params[1] + 50000))

        vel_mask_co = np.logical_and(co_specaxis.value >=
                                     max(np.min(hi_specaxis.value), hi_params[1] - 50000),
                                     co_specaxis.value <=
                                     min(np.max(hi_specaxis.value), hi_params[1] + 50000))

        # Plot to see how it's doing
        ax1 = ax[0]
        ax2 = ax1.twinx()
        ln1 = ax2.plot(hi_specaxis[vel_mask] / 1000., co_model(hi_specaxis[vel_mask]), label='CO model',
                       color=cpal[1])
        ln2 = ax1.plot(hi_specaxis[vel_mask] / 1000., hi_model(hi_specaxis[vel_mask]), label='HI model',
                       color=cpal[2])
        ax1.grid()
        ax1.set_ylabel(r"$T_{\rm HI}$ (K)")
        ax2.set_ylabel(r"$T_{\rm CO}$ (K)")
        lns = ln1 + ln2
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc='upper left', frameon=True)
        # pl.legend(frameon=True)
        # pl.xlim([np.min(hi_specaxis.value), np.max(hi_specaxis.value)])
        ax1.set_xlim([max(np.min(hi_specaxis.value) / 1000.,
                          (hi_params[1] - 50000) / 1000.),
                      min(np.max(hi_specaxis.value) / 1000.,
                          (hi_params[1] + 50000) / 1000.)])

        ax3 = ax[1]
        ax3.plot(co_specaxis[vel_mask_co] / 1000., co_spectrum.value[vel_mask_co], drawstyle='steps-mid')
        ax3.plot(hi_specaxis[vel_mask] / 1000., co_model(hi_specaxis[vel_mask]), label='CO model')
        ax3.axvline(co_params[1] - co_hwhm, color=cpal[1], linestyle='--')
        ax3.axvline(co_params[1] + co_hwhm, color=cpal[1], linestyle='--')
        ax3.set_ylabel(r"$T_{\rm CO}$ (K)")
        ax3.grid()

        ax4 = ax[2]
        ax4.plot(hi_specaxis[vel_mask] / 1000., hi_spectrum[vel_mask].value, drawstyle='steps-mid')
        ax4.plot(hi_specaxis[vel_mask] / 1000., hi_model(hi_specaxis[vel_mask]),
                 color=cpal[2])
        ln_nm = ax4.plot(hi_specaxis[vel_mask] / 1000., hi_model_nomask(hi_specaxis[vel_mask]),
                         label='HI model\n(no mask)', color=cpal[3],
                         linewidth=3, alpha=0.8, linestyle=':')
        ax4.set_ylabel(r"$T_{\rm HI}$ (K)")
        ax4.legend(loc='upper left', frameon=True)

        ax4.axvline(hi_hwhm[0], color=cpal[2], linestyle='--')
        ax4.axvline(hi_hwhm[1], color=cpal[2], linestyle='--')
        ax4.grid()
        ax4.set_xlabel("Velocity (km/s)")

        pl.tight_layout()
        pl.subplots_adjust(hspace=0.0)

        pl.draw()

        if plot_spectra:
            print([(val, stderr) for val, stderr in zip(hi_params, hi_stderrs)])
            print([(val, stderr) for val, stderr in zip(hi_params_nomask,
                                                        hi_stderrs_nomask)])
            print([(val, stderr) for val, stderr in zip(co_params, co_stderrs)])

            raw_input("{0}, {1}: {2}, {3}. Flagged: {5}: {4}"
                      .format(y, x, round(co_params[-1]),
                              round(hi_params[-1]),
                              results['multicomp_flag_HI'][i],
                              results['multicomp_flag_CO'][i]))
        else:
            # Save the examples
            fig.savefig(osjoin(figure_folder, "co_hi_spectra_models_{0}_{1}.png").format(y, x),
                        overwrite=True)
            fig.savefig(osjoin(figure_folder, "co_hi_spectra_models_{0}_{1}.pdf").format(y, x),
                        overwrite=True)

            # Record their parameters to check
            example_params.append([hi_params, hi_params_nomask, co_params])

        pl.close()


# Make a radius array
radii = gal.radius(header=cube.header).to(u.kpc)
radii_pts = radii[yposns, xposns]

# And the position angles
pang = gal.position_angles(header=cube.header).to(u.deg)
pang_pts = pang[yposns, xposns]

skycoord_grid = gal.skycoord_grid(header=cube.header)
skycoord_pts = skycoord_grid[yposns, xposns]
ra_pts = skycoord_pts.ra
dec_pts = skycoord_pts.dec

results["Rgal"] = radii_pts
results["PAgal"] = pang_pts
results["RA"] = ra_pts
results["Dec"] = dec_pts
results["ypts"] = yposns
results["xpts"] = xposns

if not os.path.exists(fourteenB_HI_data_wGBT_path("tables", no_check=True)):
    os.mkdir(fourteenB_HI_data_wGBT_path("tables", no_check=True))

# Save the lists of points in a table
tab = Table([results[key] for key in results],
            names=results.keys())

tab.write(fourteenB_HI_data_wGBT_path("tables/hi_co_gaussfit_column_densities_perpix.fits",
                                      no_check=True),
          overwrite=True)
