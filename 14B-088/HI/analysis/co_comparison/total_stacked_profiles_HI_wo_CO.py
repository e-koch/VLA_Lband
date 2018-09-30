
'''
The avg HI line width towards CO detections is larger than the stacked
line width of all HI LOS.

What are the stacked profile properties of HI where no CO is detected?
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


hi_stackpath = lambda x: osjoin(fourteenB_HI_data_wGBT_path("", no_check=True),
                                "stacked_spectra", x)
if not os.path.exists(hi_stackpath("")):
    os.mkdir(hi_stackpath(""))

dr = 500 * u.pc
max_radius = (7.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)
maxrad_string = "{0}{1}".format(int(max_radius.value), max_radius.unit)

pa_bounds_n = Angle([0.5 * np.pi * u.rad, -0.5 * np.pi * u.rad])

pa_bounds_s = Angle([-0.5 * np.pi * u.rad, 0.5 * np.pi * u.rad])

# Avg rms noise in smoothed cube is 16 mK
sigma = 2.7 * u.K

hi_cube_peakvel = \
    SpectralCube.read(fourteenB_wGBT_HI_file_dict['PeakSub_Cube'])
hi_mask = fits.open(fourteenB_wGBT_HI_file_dict['PeakSub_Mask'])[0].data > 0
hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask)

del hi_mask


co_mask = fits.open(iram_co21_14B088_data_path(
    "m33.co21_iram.14B-088_HI_source_mask.fits"))[0].data

# Require there be 3 pixels in the mask
no_co_mask_spatial = co_mask.sum(0) <= 3

# Apply this mask to the HI

hi_cube_peakvel = hi_cube_peakvel.with_mask(no_co_mask_spatial)

bin_centers, total_spectrum_hi_radial_peakvel, num_pixels = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='slice')

total_spectrum_hi_radial_peakvel_n, num_pixels_n = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='slice')[1:]

total_spectrum_hi_radial_peakvel_s, num_pixels_s = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='slice')[1:]

# Separately save the number of pixels in each bin
np.save(hi_stackpath("radial_stacking_pixelsinbin_{0}_noCO.npy").format(wstring), num_pixels)
np.save(hi_stackpath("radial_stacking_pixelsinbin_north_{0}_noCO.npy").format(wstring), num_pixels_n)
np.save(hi_stackpath("radial_stacking_pixelsinbin_south_{0}_noCO.npy").format(wstring), num_pixels_s)

spec_shape = hi_cube_peakvel.shape[0]

peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                             wcs=hi_cube_peakvel.wcs)
peakvel_stack.write(hi_stackpath("peakvel_stacked_radial_{0}_noCO.fits".format(wstring)),
                    overwrite=True)

peakvel_stack_n = SpectralCube(data=total_spectrum_hi_radial_peakvel_n.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack_n.write(hi_stackpath("peakvel_stacked_radial_north_{0}_noCO.fits".format(wstring)),
                      overwrite=True)
peakvel_stack_s = SpectralCube(data=total_spectrum_hi_radial_peakvel_s.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack_s.write(hi_stackpath("peakvel_stacked_radial_south_{0}_noCO.fits".format(wstring)),
                      overwrite=True)

total_spectrum_hi_peakvel = total_spectrum_hi_radial_peakvel.sum(0)

# Save each of these
oned_wcs = hi_cube_peakvel[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_peakvel.value,
             unit=total_spectrum_hi_peakvel.unit,
             wcs=oned_wcs).write(hi_stackpath("peakvel_stacked_{0}_noCO.fits".format(maxrad_string)),
                                 overwrite=True)

total_spectrum_hi_peakvel_n = total_spectrum_hi_radial_peakvel_n.sum(0)

# Save each of these
oned_wcs = hi_cube_peakvel[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_peakvel_n.value,
             unit=total_spectrum_hi_peakvel.unit,
             wcs=oned_wcs).write(hi_stackpath("peakvel_stacked_north_{0}_noCO.fits".format(maxrad_string)),
                                 overwrite=True)


total_spectrum_hi_peakvel_s = total_spectrum_hi_radial_peakvel_s.sum(0)

# Save each of these
oned_wcs = hi_cube_peakvel[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_peakvel_s.value,
             unit=total_spectrum_hi_peakvel.unit,
             wcs=oned_wcs).write(hi_stackpath("peakvel_stacked_south_{0}_noCO.fits".format(maxrad_string)),
                                 overwrite=True)

del hi_cube_peakvel

# Check how properties change with the masking
sigma_noise = (16 * u.mK).to(u.K).value
npix_beam = 41.


filename = hi_stackpath("peakvel_stacked_{0}_noCO.fits".format(maxrad_string))
filename_n = hi_stackpath("peakvel_stacked_north_{0}_noCO.fits".format(maxrad_string))
filename_s = hi_stackpath("peakvel_stacked_south_{0}_noCO.fits".format(maxrad_string))

stacks = OneDSpectrum.from_hdu(fits.open(filename))
stacks_n = OneDSpectrum.from_hdu(fits.open(filename_n))
stacks_s = OneDSpectrum.from_hdu(fits.open(filename_s))

filename = hi_stackpath("radial_stacking_pixelsinbin_{0}_noCO.npy".format(wstring))
filename_n = hi_stackpath("radial_stacking_pixelsinbin_north_{0}_noCO.npy".format(wstring))
filename_s = hi_stackpath("radial_stacking_pixelsinbin_south_{0}_noCO.npy".format(wstring))

num_pix = np.load(filename).sum()
num_pix_n = np.load(filename_n).sum()
num_pix_s = np.load(filename_s).sum()

# Limit the velocity range the model is fit to. The no-mask and 1-sigma
# cases have some low level systematic variations at large velocities
vels = stacks.spectral_axis.value / 1000.
vel_mask = np.logical_and(vels < 50, vels > -50)

props_all = fit_hwhm(stacks.spectral_axis[vel_mask].value,
                     stacks[vel_mask].value, asymm='full',
                     sigma_noise=sigma_noise,
                     nbeams=num_pix / npix_beam,
                     niters=1000)
props_n = fit_hwhm(stacks_n.spectral_axis[vel_mask].value,
                   stacks_n[vel_mask].value, asymm='full',
                   sigma_noise=sigma_noise,
                   nbeams=num_pix_n / npix_beam,
                   niters=1000)
props_s = fit_hwhm(stacks_s.spectral_axis[vel_mask].value,
                   stacks_s[vel_mask].value, asymm='full',
                   sigma_noise=sigma_noise,
                   nbeams=num_pix_s / npix_beam,
                   niters=1000)

# Plot shape parameters as a function of sigma

# levels = [0, 1, 2, 3]
# y_lim = [-1, 4]

# param_names = [r'$\sigma$ (m/s)', r'$v_0$ (m/s)', r'$f_{\rm wings}$',
#                r'$\sigma_{\rm wings}$ (m/s)', r'$a$', r'$\kappa$']

# twocolumn_twopanel_figure()

# fig, axs = plt.subplots(2, 3, sharex=True, figsize=(8.72, 4.75))

# for i, ax in enumerate(axs.ravel()):

#     param_vals = [props_all[j][0][i] for j in range(4)]
#     param_lows = [props_all[j][1][0][i] for j in range(4)]
#     param_highs = [props_all[j][1][1][i] for j in range(4)]

#     ax.errorbar(levels, param_vals, yerr=[param_lows, param_highs],
#                 label='All')
#     ax.set_ylabel(param_names[i])
#     ax.grid()

#     param_vals = [props_n[j][0][i] for j in range(4)]
#     param_lows = [props_n[j][1][0][i] for j in range(4)]
#     param_highs = [props_n[j][1][1][i] for j in range(4)]

#     ax.errorbar(levels, param_vals, yerr=[param_lows, param_highs],
#                 label='North')

#     param_vals = [props_s[j][0][i] for j in range(4)]
#     param_lows = [props_s[j][1][0][i] for j in range(4)]
#     param_highs = [props_s[j][1][1][i] for j in range(4)]

#     ax.errorbar(levels, param_vals, yerr=[param_lows, param_highs],
#                 label='South')

#     if i == 4:
#         ax.set_xlabel('N-sigma mask cut')

#         ax.legend(frameon=True)

# plt.tight_layout()

# fig.savefig(osjoin(figure_folder, "CO21_peakvel_HWHM_vs_masking.pdf"))
# fig.savefig(osjoin(figure_folder, "CO21_peakvel_HWHM_vs_masking.png"))

# plt.close()

default_figure()
