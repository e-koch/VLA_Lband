
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
from pandas import DataFrame
import pandas as pd
from astropy.utils.console import ProgressBar

from cube_analysis.spectral_stacking import radial_stacking
from cube_analysis.spectral_stacking_models import fit_hwhm

from paths import (fourteenB_wGBT_HI_file_dict,
                   iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path,
                   fourteenB_HI_data_path,
                   allfigs_path)
from constants import hi_freq
from galaxy_params import gal_feath as gal
from plotting_styles import default_figure, twocolumn_twopanel_figure, onecolumn_figure


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

nbins = np.int(np.floor(max_radius / dr))
inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)

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
no_co_mask_spatial = co_mask.sum(0) < 1

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
peakvel_stack = peakvel_stack.with_mask(np.ones_like(peakvel_stack, dtype='bool'))
peakvel_stack.write(hi_stackpath("peakvel_stacked_radial_{0}_noCO.fits".format(wstring)),
                    overwrite=True)

peakvel_stack_n = SpectralCube(data=total_spectrum_hi_radial_peakvel_n.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack_n = peakvel_stack_n.with_mask(np.ones_like(peakvel_stack_n, dtype='bool'))
peakvel_stack_n.write(hi_stackpath("peakvel_stacked_radial_north_{0}_noCO.fits".format(wstring)),
                      overwrite=True)
peakvel_stack_s = SpectralCube(data=total_spectrum_hi_radial_peakvel_s.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack_s = peakvel_stack_s.with_mask(np.ones_like(peakvel_stack_s, dtype='bool'))
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
sigma_noise = 2.7 * u.K
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

hi_mod_vals = {}

hi_mod_vals["fulldisk"] = props_all[0]
hi_mod_vals["fulldisk_lower"] = props_all[1][0]
hi_mod_vals["fulldisk_upper"] = props_all[1][1]

hi_mod_vals["north"] = props_n[0]
hi_mod_vals["north_lower"] = props_n[1][0]
hi_mod_vals["north_upper"] = props_n[1][1]

hi_mod_vals["south"] = props_s[0]
hi_mod_vals["south_lower"] = props_s[1][0]
hi_mod_vals["south_upper"] = props_s[1][1]

# Save these stacking test parameters in a nice table

parnames_hwhm = props_all[2]

df_hi = DataFrame(hi_mod_vals, index=parnames_hwhm)
df_hi.to_csv(fourteenB_HI_data_wGBT_path("tables/hi_totalprof_hwhm_noCO_stacking.csv",
                                         no_check=True))

# Essentially find that the line width is the same w/ or w/o LOS w/ CO
# detections. The suggests a large fraction of LOS have narrower line widths
# than those towards CO detections.

# But CO detections have a strong radial gradient. Is that part of what we're
# seeing?

param_names = parnames_hwhm

hi_radial_params = {}

file_labels = ["fulldisk", "north", "south"]

for sub in file_labels:
    for name in param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_lowlim = "{}_low_lim".format(par_name)
        par_uplim = "{}_up_lim".format(par_name)

        hi_radial_params[par_name] = np.zeros_like(inneredge.value)
        hi_radial_params[par_lowlim] = np.zeros_like(inneredge.value)
        hi_radial_params[par_uplim] = np.zeros_like(inneredge.value)


for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):
    print("On {0} of {1}".format(ctr + 1, len(inneredge)))
    hi_spectra = [peakvel_stack[:, ctr, 0],
                  peakvel_stack_n[:, ctr, 0],
                  peakvel_stack_s[:, ctr, 0]]

    for spectrum, label in ProgressBar(zip(hi_spectra, file_labels)):

        vels = spectrum.spectral_axis.to(u.km / u.s).value

        if "north" in label:
            nbeams = num_pixels_n[ctr] / npix_beam
        elif "south" in label:
            nbeams = num_pixels_s[ctr] / npix_beam
        else:
            nbeams = num_pixels[ctr] / npix_beam

        # Fit +/- 60 km/s
        vel_mask = np.logical_and(vels >= -60, vels <= 60)

        parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_HI_hwhm = \
            fit_hwhm(vels[vel_mask], spectrum.value[vel_mask],
                     sigma_noise=sigma_noise,
                     nbeams=nbeams, niters=100, interp_factor=1.)[:-1]

        for idx, name in enumerate(parnames_hwhm):
            par_name = "{0}_{1}".format(label, name)
            hi_radial_params[par_name][ctr] = parvals_hwhm[idx]
            hi_radial_params["{}_low_lim".format(par_name)][ctr] = \
                np.abs(parerrs_hwhm[0, idx])
            hi_radial_params["{}_up_lim".format(par_name)][ctr] = \
                np.abs(parerrs_hwhm[1, idx])

bin_names = ["{}-{}".format(r0.value, r1)
             for r0, r1 in zip(inneredge, outeredge)]

bin_cents = (outeredge - dr / 2.).to(u.kpc).value
hi_radial_params["bin_center"] = bin_cents

hi_radial_fits = DataFrame(hi_radial_params, index=bin_names)

# hi_radial_fits.to_latex(alltables_path("hi_gaussian_hwhm_totalprof_fits_radial.tex"))
hi_radial_fits.to_csv(fourteenB_HI_data_wGBT_path("tables/hi_gaussian_hwhm_totalprof_fits_noCO_radial.csv",
                                             no_check=True))


# Load in the non-masked fit values
hi_radial_fits_all = pd.read_csv(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_radial_500pc.csv"),
                                 index_col=0)

# Plot line width as function of radius

onecolumn_figure()

plt.errorbar(hi_radial_fits['bin_center'], hi_radial_fits_all['peaksub_sigma'],
             yerr=hi_radial_fits_all['peaksub_sigma_up_lim'], drawstyle='steps-mid',
             label='All LOS')
plt.errorbar(hi_radial_fits['bin_center'], hi_radial_fits['fulldisk_sigma'],
             yerr=hi_radial_fits['fulldisk_sigma_up_lim'], drawstyle='steps-mid',
             label='No CO LOS')
plt.legend(frameon=True)
plt.grid()
plt.ylabel(r"$\sigma_{\rm HWHM}$ (km/s)")
plt.xlabel("Radius (kpc)")
plt.ylim([4, 10])
plt.tight_layout()

# No. The no CO line width is marginally larger (~200 m/s) at some radii, but
# not vastly different than the all LOS line widths.

plt.savefig(osjoin(figure_folder, "HI_sigma_stacked_noCO_comparison.pdf"))
plt.savefig(osjoin(figure_folder, "HI_sigma_stacked_noCO_comparison.png"))

plt.close()

plt.errorbar(hi_radial_fits['bin_center'],
             hi_radial_fits_all['peaksub_f_wings'],
             yerr=hi_radial_fits_all['peaksub_f_wings_up_lim'],
             drawstyle='steps-mid',
             label='All LOS')
plt.errorbar(hi_radial_fits['bin_center'], hi_radial_fits['fulldisk_f_wings'],
             yerr=hi_radial_fits['fulldisk_f_wings_up_lim'],
             drawstyle='steps-mid',
             label='No CO LOS')
plt.legend(frameon=True)
plt.grid()
plt.ylabel(r"$f_{\rm wings}$")
plt.xlabel("Radius (kpc)")
plt.ylim([0., 0.5])
plt.tight_layout()

plt.savefig(osjoin(figure_folder, "HI_fwings_stacked_noCO_comparison.pdf"))
plt.savefig(osjoin(figure_folder, "HI_fwings_stacked_noCO_comparison.png"))

plt.close()

# How about the CO stacked w/o the bright LOS?

co_stackpath = lambda x: osjoin(iram_co21_14B088_data_path("", no_check=True),
                                "stacked_spectra", x)
if not os.path.exists(co_stackpath("")):
    os.mkdir(co_stackpath(""))

co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_feather.peakvels_corrected.fits"))

# Cut out most LOS with measurable CO emission
no_co_mask_spatial = co_mask.sum(0) < 1
co_cube = co_cube.with_mask(no_co_mask_spatial)

sigma_noise = (16 * u.mK).to(u.K).value


bin_centers, total_spectrum_co_radial, num_pixels = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='slice')

bin_centers, total_spectrum_co_radial_n, num_pixels_n = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='slice')

bin_centers, total_spectrum_co_radial_s, num_pixels_s = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='slice')

# Separately save the number of pixels in each bin
np.save(co_stackpath("radial_stacking_pixelsinbin_{0}_noCO.npy").format(wstring), num_pixels)
np.save(co_stackpath("radial_stacking_pixelsinbin_north_{0}_noCO.npy").format(wstring), num_pixels_n)
np.save(co_stackpath("radial_stacking_pixelsinbin_south_{0}_noCO.npy").format(wstring), num_pixels_s)

spec_shape = co_cube.shape[0]

peakvel_stack = SpectralCube(data=total_spectrum_co_radial.T.reshape((spec_shape, bin_centers.size, 1)),
                             wcs=co_cube.wcs)
peakvel_stack = peakvel_stack.with_mask(np.ones_like(peakvel_stack, dtype='bool'))
peakvel_stack.write(co_stackpath("peakvel_stacked_radial_{0}_noCO.fits".format(wstring)),
                    overwrite=True)

peakvel_stack_n = SpectralCube(data=total_spectrum_co_radial_n.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=co_cube.wcs)
peakvel_stack_n = peakvel_stack_n.with_mask(np.ones_like(peakvel_stack_n, dtype='bool'))
peakvel_stack_n.write(co_stackpath("peakvel_stacked_radial_north_{0}_noCO.fits".format(wstring)),
                      overwrite=True)
peakvel_stack_s = SpectralCube(data=total_spectrum_co_radial_s.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=co_cube.wcs)
peakvel_stack_s = peakvel_stack_s.with_mask(np.ones_like(peakvel_stack_s, dtype='bool'))
peakvel_stack_s.write(co_stackpath("peakvel_stacked_radial_south_{0}_noCO.fits".format(wstring)),
                      overwrite=True)

total_spectrum_co_peakvel = total_spectrum_co_radial.sum(0)

# Save each of these
oned_wcs = co_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co_peakvel.value,
             unit=total_spectrum_co_peakvel.unit,
             wcs=oned_wcs).write(co_stackpath("peakvel_stacked_{0}_noCO.fits".format(maxrad_string)),
                                 overwrite=True)

total_spectrum_co_peakvel_n = total_spectrum_co_radial_n.sum(0)

# Save each of these
oned_wcs = co_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co_peakvel_n.value,
             unit=total_spectrum_co_peakvel.unit,
             wcs=oned_wcs).write(co_stackpath("peakvel_stacked_north_{0}_noCO.fits".format(maxrad_string)),
                                 overwrite=True)


total_spectrum_co_peakvel_s = total_spectrum_co_radial_s.sum(0)

# Save each of these
oned_wcs = co_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co_peakvel_s.value,
             unit=total_spectrum_co_peakvel.unit,
             wcs=oned_wcs).write(co_stackpath("peakvel_stacked_south_{0}_noCO.fits".format(maxrad_string)),
                                 overwrite=True)

del co_cube

npix_beam = 41.


filename = co_stackpath("peakvel_stacked_{0}_noCO.fits".format(maxrad_string))
filename_n = co_stackpath("peakvel_stacked_north_{0}_noCO.fits".format(maxrad_string))
filename_s = co_stackpath("peakvel_stacked_south_{0}_noCO.fits".format(maxrad_string))

stacks = OneDSpectrum.from_hdu(fits.open(filename))
stacks_n = OneDSpectrum.from_hdu(fits.open(filename_n))
stacks_s = OneDSpectrum.from_hdu(fits.open(filename_s))

filename = co_stackpath("radial_stacking_pixelsinbin_{0}_noCO.npy".format(wstring))
filename_n = co_stackpath("radial_stacking_pixelsinbin_north_{0}_noCO.npy".format(wstring))
filename_s = co_stackpath("radial_stacking_pixelsinbin_south_{0}_noCO.npy".format(wstring))

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

co_mod_vals = {}

co_mod_vals["fulldisk"] = props_all[0]
co_mod_vals["fulldisk_lower"] = props_all[1][0]
co_mod_vals["fulldisk_upper"] = props_all[1][1]

co_mod_vals["north"] = props_n[0]
co_mod_vals["north_lower"] = props_n[1][0]
co_mod_vals["north_upper"] = props_n[1][1]

co_mod_vals["south"] = props_s[0]
co_mod_vals["south_lower"] = props_s[1][0]
co_mod_vals["south_upper"] = props_s[1][1]

# Save these stacking test parameters in a nice table

parnames_hwhm = props_all[2]

df_co = DataFrame(co_mod_vals, index=parnames_hwhm)
df_co.to_csv(iram_co21_14B088_data_path("tables/co_totalprof_hwhm_noCO_stacking.csv",
                                        no_check=True))

# Essentially find that the line width is the same w/ or w/o LOS w/ CO
# detections. The suggests a large fraction of LOS have narrower line widths
# than those towards CO detections.

# But CO detections have a strong radial gradient. Is that part of what we're
# seeing?

param_names = parnames_hwhm

co_radial_params = {}

file_labels = ["fulldisk", "north", "south"]

for sub in file_labels:
    for name in param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_lowlim = "{}_low_lim".format(par_name)
        par_uplim = "{}_up_lim".format(par_name)

        co_radial_params[par_name] = np.zeros_like(inneredge.value)
        co_radial_params[par_lowlim] = np.zeros_like(inneredge.value)
        co_radial_params[par_uplim] = np.zeros_like(inneredge.value)


for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):
    print("On {0} of {1}".format(ctr + 1, len(inneredge)))
    co_spectra = [peakvel_stack[:, ctr, 0],
                  peakvel_stack_n[:, ctr, 0],
                  peakvel_stack_s[:, ctr, 0]]

    for spectrum, label in ProgressBar(zip(co_spectra, file_labels)):

        vels = spectrum.spectral_axis.to(u.km / u.s).value

        if "north" in label:
            nbeams = num_pixels_n[ctr] / npix_beam
        elif "south" in label:
            nbeams = num_pixels_s[ctr] / npix_beam
        else:
            nbeams = num_pixels[ctr] / npix_beam

        # Fit +/- 60 km/s
        vel_mask = np.logical_and(vels >= -60, vels <= 60)

        parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_CO_hwhm = \
            fit_hwhm(vels[vel_mask], spectrum.value[vel_mask],
                     sigma_noise=sigma_noise,
                     nbeams=nbeams, niters=100, interp_factor=1.)[:-1]

        for idx, name in enumerate(parnames_hwhm):
            par_name = "{0}_{1}".format(label, name)
            co_radial_params[par_name][ctr] = parvals_hwhm[idx]
            co_radial_params["{}_low_lim".format(par_name)][ctr] = \
                np.abs(parerrs_hwhm[0, idx])
            co_radial_params["{}_up_lim".format(par_name)][ctr] = \
                np.abs(parerrs_hwhm[1, idx])

bin_names = ["{}-{}".format(r0.value, r1)
             for r0, r1 in zip(inneredge, outeredge)]

bin_cents = (outeredge - dr / 2.).to(u.kpc).value
co_radial_params["bin_center"] = bin_cents

co_radial_fits = DataFrame(co_radial_params, index=bin_names)

co_radial_fits.to_csv(iram_co21_14B088_data_path("tables/co_gaussian_hwhm_totalprof_fits_noCO_radial.csv",
                                                  no_check=True))


# Load in the non-masked fit values
co_radial_fits_all = pd.read_csv(iram_co21_14B088_data_path("tables/co_hwhm_totalprof_fits_radial_500pc.csv"),
                                 index_col=0)

# Plot line width as function of radius

onecolumn_figure()

plt.errorbar(co_radial_fits['bin_center'], co_radial_fits_all['peaksub_sigma'],
             yerr=co_radial_fits_all['peaksub_sigma_up_lim'], drawstyle='steps-mid',
             label='All LOS')
plt.errorbar(co_radial_fits['bin_center'], co_radial_fits['fulldisk_sigma'],
             yerr=co_radial_fits['fulldisk_sigma_up_lim'], drawstyle='steps-mid',
             label='No CO LOS')
plt.legend(frameon=True)
plt.grid()
plt.ylabel(r"$\sigma_{\rm HWHM}$ (km/s)")
plt.xlabel("Radius (kpc)")
plt.tight_layout()

# No. The no CO line width is marginally larger (~200 m/s) at some radii, but
# not vastly different than the all LOS line widths.

plt.savefig(osjoin(figure_folder, "CO_sigma_stacked_noCO_comparison.pdf"))
plt.savefig(osjoin(figure_folder, "CO_sigma_stacked_noCO_comparison.png"))

plt.close()

plt.errorbar(co_radial_fits['bin_center'],
             co_radial_fits_all['peaksub_f_wings'],
             yerr=co_radial_fits_all['peaksub_f_wings_up_lim'],
             drawstyle='steps-mid',
             label='All LOS')
plt.errorbar(co_radial_fits['bin_center'], co_radial_fits['fulldisk_f_wings'],
             yerr=co_radial_fits['fulldisk_f_wings_up_lim'],
             drawstyle='steps-mid',
             label='No CO LOS')
plt.legend(frameon=True)
plt.grid()
plt.ylabel(r"$f_{\rm wings}$")
plt.xlabel("Radius (kpc)")
plt.ylim([-0.25, 0.5])
plt.tight_layout()

plt.savefig(osjoin(figure_folder, "CO_fwings_stacked_noCO_comparison.pdf"))
plt.savefig(osjoin(figure_folder, "CO_fwings_stacked_noCO_comparison.png"))

plt.close()

default_figure()
