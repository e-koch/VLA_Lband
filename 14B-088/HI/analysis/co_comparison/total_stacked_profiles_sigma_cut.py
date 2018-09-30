
'''
Repeat the peak stacking from total_stacked_profiles, but impose a brightness
cut on the CO.

This is to:

1) determine if the line wing strength scales with number of LOS, which
should be evident if it arises from the error beam.

2) investigate if the asymmetry persists in a limited number of LOS.

3) If the key properties (e.g. line width) change with a different sample
(unlikely based on the stack of spectra used in the LOS analysis).

4) Determine importance of the high S/N peaks that were removed from the
single Gaussian good spectra.
'''

import astropy.units as u
from spectral_cube import SpectralCube, OneDSpectrum
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from os.path import join as osjoin
import os
from astropy.coordinates import Angle

from cube_analysis.spectral_stacking import radial_stacking
from cube_analysis.spectral_stacking_models import fit_hwhm

from paths import (fourteenB_wGBT_HI_file_dict,
                   iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path,
                   allfigs_path)
from constants import hi_freq
from galaxy_params import gal_feath as gal
from plotting_styles import default_figure, twocolumn_twopanel_figure


default_figure()

figure_folder = allfigs_path("stacked_profiles")
if not os.path.exists(figure_folder):
    os.mkdir(figure_folder)


co_stackpath = lambda x: osjoin(iram_co21_14B088_data_path("", no_check=True), "stacked_spectra", x)
if not os.path.exists(co_stackpath("")):
    os.mkdir(co_stackpath(""))

dr = 500 * u.pc
max_radius = (7.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)
maxrad_string = "{0}{1}".format(int(max_radius.value), max_radius.unit)

pa_bounds_n = Angle([0.5 * np.pi * u.rad, -0.5 * np.pi * u.rad])

pa_bounds_s = Angle([-0.5 * np.pi * u.rad, 0.5 * np.pi * u.rad])

# Avg rms noise in smoothed cube is 16 mK
sigma = 16. * u.mK

co_cube_peakvel = \
    SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_feather.peakvels_corrected.fits"))
co_cube_peakvel.allow_huge_operations = True

for level in range(1, 4):

    mask = co_cube_peakvel > level * sigma

    # Require there be >= 3 pixels above the threshold to be included
    spat_mask = mask.include().sum(0) >= 3

    masked_cube = co_cube_peakvel.with_mask(spat_mask)

    bin_centers, total_spectrum_co_radial_peakvel, num_pixels = \
        radial_stacking(gal, masked_cube, dr=dr,
                        max_radius=max_radius,
                        pa_bounds=None,
                        verbose=True,
                        how='cube')

    total_spectrum_co_radial_peakvel_n, num_pixels_n = \
        radial_stacking(gal, masked_cube, dr=dr,
                        max_radius=max_radius,
                        pa_bounds=pa_bounds_n,
                        verbose=True,
                        how='cube')[1:]

    total_spectrum_co_radial_peakvel_s, num_pixels_s = \
        radial_stacking(gal, masked_cube, dr=dr,
                        max_radius=max_radius,
                        pa_bounds=pa_bounds_s,
                        verbose=True,
                        how='cube')[1:]

    # Separately save the number of pixels in each bin
    np.save(co_stackpath("radial_stacking_pixelsinbin_{0}_sigmacut_{1}.npy").format(wstring, level), num_pixels)
    np.save(co_stackpath("radial_stacking_pixelsinbin_north_{0}_sigmacut_{1}.npy").format(wstring, level), num_pixels_n)
    np.save(co_stackpath("radial_stacking_pixelsinbin_south_{0}_sigmacut_{1}.npy").format(wstring, level), num_pixels_s)

    spec_shape = co_cube_peakvel.shape[0]

    peakvel_stack = SpectralCube(data=total_spectrum_co_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                                 wcs=co_cube_peakvel.wcs)
    peakvel_stack.write(co_stackpath("peakvel_stacked_radial_{0}_sigmacut_{1}.fits".format(wstring, level)),
                        overwrite=True)

    peakvel_stack_n = SpectralCube(data=total_spectrum_co_radial_peakvel_n.T.reshape((spec_shape, bin_centers.size, 1)),
                                   wcs=co_cube_peakvel.wcs)
    peakvel_stack_n.write(co_stackpath("peakvel_stacked_radial_north_{0}_sigmacut_{1}.fits".format(wstring, level)),
                          overwrite=True)
    peakvel_stack_s = SpectralCube(data=total_spectrum_co_radial_peakvel_s.T.reshape((spec_shape, bin_centers.size, 1)),
                                   wcs=co_cube_peakvel.wcs)
    peakvel_stack_s.write(co_stackpath("peakvel_stacked_radial_south_{0}_sigmacut_{1}.fits".format(wstring, level)),
                          overwrite=True)

    total_spectrum_co_peakvel = total_spectrum_co_radial_peakvel.sum(0)

    # Save each of these
    oned_wcs = co_cube_peakvel[:, 0, 0].wcs
    OneDSpectrum(total_spectrum_co_peakvel.value,
                 unit=total_spectrum_co_peakvel.unit,
                 wcs=oned_wcs).write(co_stackpath("peakvel_stacked_{0}_sigmacut_{1}.fits".format(maxrad_string, level)),
                                     overwrite=True)

    total_spectrum_co_peakvel_n = total_spectrum_co_radial_peakvel_n.sum(0)

    # Save each of these
    oned_wcs = co_cube_peakvel[:, 0, 0].wcs
    OneDSpectrum(total_spectrum_co_peakvel_n.value,
                 unit=total_spectrum_co_peakvel.unit,
                 wcs=oned_wcs).write(co_stackpath("peakvel_stacked_north_{0}_sigmacut_{1}.fits".format(maxrad_string, level)),
                                     overwrite=True)


    total_spectrum_co_peakvel_s = total_spectrum_co_radial_peakvel_s.sum(0)

    # Save each of these
    oned_wcs = co_cube_peakvel[:, 0, 0].wcs
    OneDSpectrum(total_spectrum_co_peakvel_s.value,
                 unit=total_spectrum_co_peakvel.unit,
                 wcs=oned_wcs).write(co_stackpath("peakvel_stacked_south_{0}_sigmacut_{1}.fits".format(maxrad_string, level)),
                                     overwrite=True)
del co_cube_peakvel

# Save version of the no masking radial stacking in the N and S
total_spectrum_co_radial = SpectralCube.read(co_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)))
total_spectrum_co_radial_n = SpectralCube.read(co_stackpath("peakvel_stacked_radial_north_{}.fits".format(wstring)))
total_spectrum_co_radial_s = SpectralCube.read(co_stackpath("peakvel_stacked_radial_south_{}.fits".format(wstring)))

total_spectrum_co_radial_n.sum(axis=(1, 2)).write(co_stackpath("peakvel_stacked_north_{0}.fits".format(maxrad_string)),
                                                  overwrite=True)
total_spectrum_co_radial_s.sum(axis=(1, 2)).write(co_stackpath("peakvel_stacked_south_{0}.fits".format(maxrad_string)),
                                                  overwrite=True)

# Check how properties change with the masking
props_all = []
props_n = []
props_s = []

stacks = []
stacks_n = []
stacks_s = []

num_pix = []
num_pix_n = []
num_pix_s = []

sigma_noise = (16 * u.mK).to(u.K).value
npix_beam = 41.


for level in range(0, 4):

    if level == 0:
        filename = co_stackpath("peakvel_stacked_{0}.fits".format(maxrad_string))
        filename_n = co_stackpath("peakvel_stacked_north_{0}.fits".format(maxrad_string))
        filename_s = co_stackpath("peakvel_stacked_south_{0}.fits".format(maxrad_string))
    else:
        filename = co_stackpath("peakvel_stacked_{0}_sigmacut_{1}.fits".format(maxrad_string, level))
        filename_n = co_stackpath("peakvel_stacked_north_{0}_sigmacut_{1}.fits".format(maxrad_string, level))
        filename_s = co_stackpath("peakvel_stacked_south_{0}_sigmacut_{1}.fits".format(maxrad_string, level))

    stacks.append(OneDSpectrum.from_hdu(fits.open(filename)))
    stacks_n.append(OneDSpectrum.from_hdu(fits.open(filename_n)))
    stacks_s.append(OneDSpectrum.from_hdu(fits.open(filename_s)))

    if level == 0:
        filename = co_stackpath("radial_stacking_pixelsinbin_{0}.npy".format(wstring))
        filename_n = co_stackpath("radial_stacking_pixelsinbin_north_{0}.npy".format(wstring))
        filename_s = co_stackpath("radial_stacking_pixelsinbin_south_{0}.npy".format(wstring))
    else:
        filename = co_stackpath("radial_stacking_pixelsinbin_{0}_sigmacut_{1}.npy".format(wstring, level))
        filename_n = co_stackpath("radial_stacking_pixelsinbin_north_{0}_sigmacut_{1}.npy".format(wstring, level))
        filename_s = co_stackpath("radial_stacking_pixelsinbin_south_{0}_sigmacut_{1}.npy".format(wstring, level))

    num_pix.append(np.load(filename).sum())
    num_pix_n.append(np.load(filename_n).sum())
    num_pix_s.append(np.load(filename_s).sum())

    # Limit the velocity range the model is fit to. The no-mask and 1-sigma
    # cases have some low level systematic variations at large velocities
    vels = stacks[level].spectral_axis.value / 1000.
    vel_mask = np.logical_and(vels < 50, vels > -50)

    props_all.append(fit_hwhm(stacks[level].spectral_axis[vel_mask].value,
                              stacks[level][vel_mask].value, asymm='full',
                              sigma_noise=sigma_noise,
                              nbeams=num_pix[level] / npix_beam,
                              niters=1000))
    props_n.append(fit_hwhm(stacks_n[level].spectral_axis[vel_mask].value,
                            stacks_n[level][vel_mask].value, asymm='full',
                            sigma_noise=sigma_noise,
                            nbeams=num_pix_n[level] / npix_beam,
                            niters=1000))
    props_s.append(fit_hwhm(stacks_s[level].spectral_axis[vel_mask].value,
                            stacks_s[level][vel_mask].value, asymm='full',
                            sigma_noise=sigma_noise,
                            nbeams=num_pix_s[level] / npix_beam,
                            niters=1000))

# Plot shape parameters as a function of sigma

levels = [0, 1, 2, 3]
y_lim = [-1, 4]

param_names = [r'$\sigma$ (m/s)', r'$v_0$ (m/s)', r'$f_{\rm wings}$',
               r'$\sigma_{\rm wings}$ (m/s)', r'$a$', r'$\kappa$']

twocolumn_twopanel_figure()

fig, axs = plt.subplots(2, 3, sharex=True, figsize=(8.72, 4.75))

for i, ax in enumerate(axs.ravel()):

    param_vals = [props_all[j][0][i] for j in range(4)]
    param_lows = [props_all[j][1][0][i] for j in range(4)]
    param_highs = [props_all[j][1][1][i] for j in range(4)]

    ax.errorbar(levels, param_vals, yerr=[param_lows, param_highs],
                label='All')
    ax.set_ylabel(param_names[i])
    ax.grid()

    param_vals = [props_n[j][0][i] for j in range(4)]
    param_lows = [props_n[j][1][0][i] for j in range(4)]
    param_highs = [props_n[j][1][1][i] for j in range(4)]

    ax.errorbar(levels, param_vals, yerr=[param_lows, param_highs],
                label='North')

    param_vals = [props_s[j][0][i] for j in range(4)]
    param_lows = [props_s[j][1][0][i] for j in range(4)]
    param_highs = [props_s[j][1][1][i] for j in range(4)]

    ax.errorbar(levels, param_vals, yerr=[param_lows, param_highs],
                label='South')

    if i == 4:
        ax.set_xlabel('N-sigma mask cut')

        ax.legend(frameon=True)

plt.tight_layout()

fig.savefig(osjoin(figure_folder, "CO21_peakvel_HWHM_vs_masking.pdf"))
fig.savefig(osjoin(figure_folder, "CO21_peakvel_HWHM_vs_masking.png"))

plt.close()

default_figure()
