
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
import os
import seaborn as sb

from cube_analysis.spectral_stacking_models import find_hwhm, fit_gaussian

from paths import (fourteenB_wGBT_HI_file_dict, iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path)
from plotting_styles import default_figure, onecolumn_figure, twocolumn_figure
from galaxy_params import gal_feath as gal
from constants import (co21_mass_conversion, hi_mass_conversion, hi_freq)


default_figure()
cpal = sb.color_palette()

# cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['PeakSub_Cube'])
cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])
co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.fits"))

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
           "coldens_HI_FWHM": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_CO_FWHM": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_HI_gauss": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "coldens_CO_gauss": np.zeros((yposns.shape)) * u.solMass / u.pc**2,
           "multicomp_flag": np.zeros(yposns.shape, dtype=bool)}

# Correct for the disk inclincation
inc = np.cos(gal.inclination)

# 30 m beam efficiency
beam_eff = 0.75

for i, (y, x) in enumerate(zip(yposns, xposns)):

    co_spectrum = co_cube[:, y, x]
    # First beam is the largest. Change is far far smaller than 1 pixel
    hi_spectrum = cube[:, y, x].to(u.K, cube.beams[0].jtok_equiv(hi_freq))

    co_params, co_cov, co_parnames, co_model = \
        fit_gaussian(co_specaxis,
                     co_spectrum.quantity)

    co_stderrs = np.sqrt(np.diag(co_cov))

    results["amp_CO"][i] = co_params[0] * u.K
    results["mean_CO"][i] = co_params[1] * u.m / u.s
    results["sigma_CO"][i] = co_params[2] * u.m / u.s
    # Adjust the mean and sigma to account for finite channel widths
    # Take to be half of the channel width.
    results["amp_stderr_CO"][i] = co_stderrs[0] * u.K
    results["mean_stderr_CO"][i] = np.sqrt(co_stderrs[1]**2 +
                                           (co_chanwidth / 2.)**2) * u.m / u.s
    results["sigma_stderr_CO"][i] = np.sqrt(co_stderrs[2]**2 +
                                            (co_chanwidth / 2.)**2) * u.m / u.s

    # Make a mask centered around the CO peak with 5x the CO FWHM
    # This will hopefully help centre the HI fit to the CO if there is
    # a brighter velocity component.
    co_hwhm = co_params[-1] * np.sqrt(2 * np.log(2))
    co_mask = np.logical_and(hi_specaxis.value >= co_params[1] - 5 * co_hwhm,
                             hi_specaxis.value <= co_params[1] + 5 * co_hwhm)

    # Limit the HI fitting to the peak. Smooth the HI spectrum to get a good
    # HWHM estimate
    # This will fail for spectra with multiple components blended together.
    # Flag those that fail and make some rough estimates instead.
    try:
        hi_sigest, hi_hwhm = find_hwhm(hi_specaxis[co_mask],
                                       hi_spectrum.spectral_smooth(kern)[co_mask])[:2]
        results["multicomp_flag"][i] = False
    except ValueError:
        # Just adopt some typical values
        hi_sigest = 7000.
        hwhm_factor = np.sqrt(2 * np.log(2))
        peak_posn = np.nanargmax(hi_spectrum.spectral_smooth(kern)[co_mask])
        hi_hwhm = np.array([hi_specaxis[co_mask][peak_posn].value - hi_sigest * hwhm_factor,
                            hi_specaxis[co_mask][peak_posn].value + hi_sigest * hwhm_factor])
        results["multicomp_flag"][i] = True

    hwhm_mask = np.logical_and(hi_specaxis.value >= hi_hwhm[0] - mask_pad,
                               hi_specaxis.value <= hi_hwhm[1] + mask_pad)

    # Give some initial guesses
    p0 = (np.nanmax(hi_spectrum), hi_specaxis[np.nanargmax(hi_spectrum)],
          hi_sigest)

    hi_params, hi_cov, hi_parnames, hi_model = \
        fit_gaussian(hi_specaxis[hwhm_mask],
                     hi_spectrum.quantity[hwhm_mask], p0=p0)

    hi_stderrs = np.sqrt(np.diag(hi_cov))

    results["amp_HI"][i] = hi_params[0] * u.K
    results["mean_HI"][i] = hi_params[1] * u.m / u.s
    results["sigma_HI"][i] = hi_params[2] * u.m / u.s
    # Adjust the mean and sigma to account for finite channel widths
    # Take to be half of the channel width.
    results["amp_stderr_HI"][i] = hi_stderrs[0] * u.K
    results["mean_stderr_HI"][i] = np.sqrt(hi_stderrs[1]**2 +
                                           (hi_chanwidth / 2.)**2) * u.m / u.s
    results["sigma_stderr_HI"][i] = np.sqrt(hi_stderrs[2]**2 +
                                            (hi_chanwidth / 2.)**2) * u.m / u.s

    # Finally calculate column densities from:
    # 1) within the FWHM of each spectra,
    # 2) from the integral of the gaussian models

    hi_hwhm_model = hi_params[-1] * np.sqrt(2 * np.log(2))
    hi_fwhm_mask = \
        np.logical_and(hi_specaxis.value >= hi_params[1] - hi_hwhm_model,
                       hi_specaxis.value <= hi_params[1] + hi_hwhm_model)
    results["coldens_HI_FWHM"][i] = \
        np.nansum(hi_spectrum[hi_fwhm_mask]) * \
        (hi_chanwidth / 1000. * u.km / u.s) * inc * hi_mass_conversion

    co_fwhm_mask = \
        np.logical_and(co_specaxis.value >= co_params[1] - co_hwhm,
                       co_specaxis.value <= co_params[1] + co_hwhm)
    results["coldens_CO_FWHM"][i] = \
        np.nansum(co_spectrum[co_fwhm_mask]) * \
        (co_chanwidth / 1000. * u.km / u.s) * \
        inc * co21_mass_conversion / beam_eff

    # 2)

    results["coldens_HI_gauss"][i] = \
        (hi_params[0] * u.K) * np.sqrt(2 * np.pi) * \
        (hi_params[-1] / 1000. * u.km / u.s) * inc * hi_mass_conversion

    results["coldens_CO_gauss"][i] = \
        (co_params[0] * u.K) * np.sqrt(2 * np.pi) * \
        (co_params[-1] / 1000. * u.km / u.s) * \
        inc * co21_mass_conversion / beam_eff

    if plot_spectra:
        # Plot to see how it's doing
        ax1 = pl.subplot(111)
        ax2 = ax1.twinx()
        ax1.plot(hi_specaxis, hi_spectrum.value, drawstyle='steps-mid')
        ax2.plot(co_specaxis, co_spectrum.value, drawstyle='steps-mid')
        ax2.plot(hi_specaxis, co_model(hi_specaxis.value), label='CO model')
        ax1.plot(hi_specaxis, hi_model(hi_specaxis.value), label='HI model')
        pl.axvline(hi_hwhm[0])
        pl.axvline(hi_hwhm[1])
        pl.legend(frameon=True)

        pl.draw()

        raw_input("{0}, {1}: {2}, {3}".format(y, x, round(co_params[-1]),
                                              round(hi_params[-1])))

        pl.clf()


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

results["Rgal"] = radii
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
