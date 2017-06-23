
import astropy.units as u
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import largest_beam
import numpy as np
import matplotlib.pyplot as p
from pandas import DataFrame
from astropy.modeling import models, fitting
from astropy.io import fits
from os.path import join
import os

from cube_analysis.spectral_stacking import total_profile

from paths import (fourteenB_HI_file_dict, iram_co21_14B088_data_path,
                   allfigs_path, alltables_path, fourteenB_HI_data_path)
from constants import co21_mass_conversion, hi_freq
from galaxy_params import gal_feath as gal
from plotting_styles import *

'''
Create profiles of HI and CO after subtracting velocity surfaces.
'''

figure_folder = "stacked_profiles"
if not os.path.exists(allfigs_path(figure_folder)):
    os.mkdir(allfigs_path(figure_folder))

dr = 500 * u.pc
max_radius = (7.0 * u.kpc).to(u.pc)

co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.rotation_corrected.fits"))
co_cube_cent = \
    SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.centroid_corrected.fits"))
co_cube_peakvel = \
    SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.peakvels_corrected.fits"))

hi_cube = SpectralCube.read(fourteenB_HI_file_dict["RotSub_Cube"])
hi_mask = fits.open(fourteenB_HI_file_dict["RotSub_Mask"])[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

hi_cube_cent = SpectralCube.read(fourteenB_HI_file_dict["CentSub_Cube"])
hi_mask_cent = fits.open(fourteenB_HI_file_dict["CentSub_Mask"])[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

hi_cube_peakvel = SpectralCube.read(fourteenB_HI_file_dict["PeakSub_Cube"])
hi_mask_peakvel = fits.open(fourteenB_HI_file_dict["PeakSub_Mask"])[0]
hi_cube_peakvel = hi_cube_cent.with_mask(hi_mask_peakvel.data > 0)

if hasattr(hi_cube, "beams"):
    hi_beam = largest_beam(hi_cube.beams)
else:
    hi_beam = hi_cube.beam

hi_radius = gal.radius(header=hi_cube.header)
co_radius = gal.radius(header=co_cube.header)

# Perform the same analysis split up into radial bins
nbins = np.int(np.floor(max_radius / dr))

inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)

total_spectrum_hi_radial = np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_co_radial = np.zeros((inneredge.size, co_cube.shape[0])) * u.K

total_spectrum_hi_radial_cent = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_co_radial_cent = \
    np.zeros((inneredge.size, co_cube.shape[0])) * u.K

total_spectrum_hi_radial_peakvel = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_co_radial_peakvel = \
    np.zeros((inneredge.size, co_cube.shape[0])) * u.K

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {}".format(r0.value, r1))

    hi_rad_mask = np.logical_and(hi_radius >= r0,
                                 hi_radius < r1)
    co_rad_mask = np.logical_and(co_radius >= r0,
                                 co_radius < r1)

    total_spectrum_hi_radial[ctr] = \
        total_profile(hi_cube, hi_rad_mask).to(u.K, hi_beam.jtok_equiv(hi_freq)).quantity

    total_spectrum_hi_radial_cent[ctr] = \
        total_profile(hi_cube_cent, hi_rad_mask).to(u.K, hi_beam.jtok_equiv(hi_freq)).quantity

    total_spectrum_hi_radial_peakvel[ctr] = \
        total_profile(hi_cube_peakvel, hi_rad_mask).to(u.K, hi_beam.jtok_equiv(hi_freq)).quantity

    total_spectrum_co_radial[ctr] = total_profile(co_cube, co_rad_mask).quantity

    total_spectrum_co_radial_cent[ctr] = \
        total_profile(co_cube_cent, co_rad_mask).quantity
    total_spectrum_co_radial_peakvel[ctr] = \
        total_profile(co_cube_peakvel, co_rad_mask).quantity

# Need to get portions of HI emission beyond 6 kpc.
total_spectrum_hi = \
    total_profile(hi_cube).to(u.K, hi_beam.jtok_equiv(hi_freq)).quantity

total_spectrum_hi_cent = \
    total_profile(hi_cube_cent).to(u.K, hi_beam.jtok_equiv(hi_freq)).quantity

total_spectrum_hi_peakvel = \
    total_profile(hi_cube_peakvel).to(u.K, hi_beam.jtok_equiv(hi_freq)).quantity

# Significant CO emission is limited to within about 6 kpc
total_spectrum_co = total_spectrum_co_radial.sum(0)

total_spectrum_co_cent = total_spectrum_co_radial_cent.sum(0)

total_spectrum_co_peakvel = total_spectrum_co_radial_peakvel.sum(0)

twocolumn_twopanel_figure()
# Plot the profiles.
fig, ax = p.subplots(1, 3, sharey=True)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)
ax[0].plot(hi_cube.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi / total_spectrum_hi.max()).value,
           'b-', drawstyle='steps-mid', label="HI")
ax[0].plot(co_cube.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co / total_spectrum_co.max()).value,
           'g--', drawstyle='steps-mid', label="CO(2-1)")
ax[0].set_xlabel("Velocity (km/s)")
ax[0].set_ylabel("Normalized Total Intensity")
# ax[0].set_title("Rotation subtracted")
ax[0].text(-45, 0.88, "Rotation\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].set_ylim([-0.02, 1.1])
ax[0].set_xlim([-50, 50])
ax[0].grid()

ax[1].plot(hi_cube_cent.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi_cent / total_spectrum_hi_cent.max()).value,
           'b-', drawstyle='steps-mid')
ax[1].plot(co_cube.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co_cent / total_spectrum_co_cent.max()).value,
           'g--', drawstyle='steps-mid')
# ax[1].set_title("Centroid subtracted")
ax[1].text(-45, 0.88, "Centroid\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[1].set_xlabel("Velocity (km/s)")
ax[1].set_ylim([-0.02, 1.1])
ax[1].set_xlim([-50, 50])
ax[1].grid()

ax[2].plot(hi_cube_peakvel.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi_peakvel / total_spectrum_hi_peakvel.max()).value,
           'b-', drawstyle='steps-mid', label="HI")
ax[2].plot(co_cube.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co_peakvel / total_spectrum_co_peakvel.max()).value,
           'g--', drawstyle='steps-mid', label="CO(2-1)")
# ax[2].set_title("Centroid subtracted")
ax[2].text(-45, 0.88, "Peak Vel.\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].set_xlabel("Velocity (km/s)")
ax[2].set_ylim([-0.02, 1.1])
ax[2].set_xlim([-50, 50])
ax[2].legend(frameon=True)
ax[2].grid()

p.tight_layout()
p.draw()

fig.savefig(allfigs_path(join(figure_folder, "total_profile_corrected_velocity_HI_CO21.pdf")))
fig.savefig(allfigs_path(join(figure_folder, "total_profile_corrected_velocity_HI_CO21.png")))

p.close()
# raw_input("Next plot?")

# Total CO mass. Using 6.7 Msol / pc^2 / K km s^-1\
pixscale = gal.distance.to(u.pc) * (np.pi / 180.) * \
    np.abs(co_cube.header["CDELT2"])
chan_width = \
    np.abs(co_cube.spectral_axis[1] - co_cube.spectral_axis[0]).to(u.km / u.s)
beam_eff = 0.75  # Beam efficiency of IRAM @ 235 GHz
# Where total_spectrum_co is in K
total_co_mass = \
    total_spectrum_co[total_spectrum_co > 0 * u.K].sum() * chan_width * \
    pixscale ** 2 * co21_mass_conversion / beam_eff
print("Total H2 Mass from CO is {} Msol".format(total_co_mass))


# We want to find the widths of these profiles.
# CO is close enough to a gaussian
# The HI has wings. Are they real wings? Maybe.. but at least some portion
# is due to double-peaked spectra. We're going to model w/ two gaussians
# with the same mean.

# HI model
g_HI_init = models.Gaussian1D(amplitude=0.75, mean=0., stddev=5.) +  \
    models.Gaussian1D(amplitude=0.25, mean=0., stddev=20.)


# Force to the same mean
def tie_mean(mod):
    return mod.mean_0


g_HI_init.mean_1.tied = tie_mean

fit_g = fitting.LevMarLSQFitter()

spectra = [total_spectrum_hi, total_spectrum_hi_cent,
           total_spectrum_hi_peakvel]
labels = ["Rot. Sub.", "Cent. Sub.", "Peak Sub."]
file_labels = ["rotsub", "centsub", "peaksub"]
hi_fit_vals = {}

vels = hi_cube.spectral_axis.to(u.km / u.s).value

onecolumn_figure(fig_ratio=0.8)

for spectrum, label, file_label in zip(spectra, labels, file_labels):

    norm_intens = (spectrum / np.nanmax(spectrum)).value

    g_HI = fit_g(g_HI_init, vels, norm_intens)

    # The covariance matrix is hidden away... tricksy
    cov = fit_g.fit_info['param_cov']
    parnames = [n for n in g_HI.param_names if n not in ['mean_1']]
    parvals = [v for (n, v) in zip(g_HI.param_names, g_HI.parameters)
               if n in parnames]

    hi_fit_vals[label + " Params"] = parvals
    hi_fit_vals[label + " Errors"] = \
        [np.sqrt(cov[i, i]) for i in range(cov.shape[0])]

    # Note that the statistical errors on the mean are too small.
    ax1 = p.subplot2grid((5, 1), (0, 0), rowspan=4)
    p.plot(vels, norm_intens, 'b-', drawstyle='steps-mid')
    p.plot(vels, g_HI(vels), 'k:', label="Total Fit")
    p.plot(vels, g_HI["None_0"](vels), 'g--', label="Narrow")
    p.plot(vels, g_HI["None_1"](vels), 'm-.', label="Wide")
    p.ylabel("HI Normalized Intensity")
    p.xticks([])
    p.xlim([-50, 50])
    p.legend(frameon=True)
    p.ylim([-0.1, 1.1])
    p.grid()

    ax2 = p.subplot2grid((5, 1), (4, 0))
    p.plot(vels, norm_intens - g_HI(vels), 'b-', drawstyle='steps-mid')
    p.hlines(0.0, -100, 100, color='k')
    p.grid()
    p.ylim([-0.015, 0.015])
    p.xlim([-50, 50])
    p.xlabel("Velocity (km/s)")
    p.ylabel("Residuals")

    p.tight_layout()
    p.subplots_adjust(hspace=0)

    filename = "total_profile_corrected_velocity_{}_hi_fit".format(file_label)
    p.savefig(allfigs_path(join(figure_folder, filename + ".pdf")))
    p.savefig(allfigs_path(join(figure_folder, filename + ".png")))
    # raw_input("Next plot?")
    p.close()

# Save parameter table
hi_param_df = DataFrame(hi_fit_vals, index=parnames)
hi_param_df.to_latex(alltables_path("hi_gaussian_totalprof_fits.tex"))
hi_param_df.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits.csv",
                                          no_check=True))

# CO model

g_CO_init = models.Gaussian1D(amplitude=1., mean=0., stddev=9.)

fit_g_co = fitting.LevMarLSQFitter()

spectra = [total_spectrum_co, total_spectrum_co_cent,
           total_spectrum_co_peakvel]
co_fit_vals = {}

co_vels = co_cube.spectral_axis.to(u.km / u.s).value

for spectrum, label, file_label in zip(spectra, labels, file_labels):

    norm_co_intens = spectrum / spectrum.max()
    g_CO = fit_g_co(g_CO_init, co_vels, norm_co_intens)

    cov = fit_g_co.fit_info['param_cov']

    co_fit_vals[label + " Params"] = g_CO.parameters
    co_fit_vals[label + " Errors"] = \
        [np.sqrt(cov[i, i]) for i in range(cov.shape[0])]

    # Better sampling for plotting
    more_vels = np.arange(co_vels.min(), co_vels.max(), 0.5)

    ax1 = p.subplot2grid((5, 1), (0, 0), rowspan=4)
    p.plot(co_vels, norm_co_intens, 'b-', drawstyle='steps-mid')
    p.plot(more_vels, g_CO(more_vels), 'k--', label="Fit")
    p.ylabel("CO Normalized Intensity")
    p.xticks([])
    p.xlim([-50, 50])
    p.legend(frameon=True)
    p.ylim([-0.1, 1.1])
    p.grid()

    ax2 = p.subplot2grid((5, 1), (4, 0))
    p.plot(co_vels, norm_co_intens - g_CO(co_vels), 'b-',
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

    filename = "total_profile_corrected_velocity_{}_co21_fit".format(file_label)
    p.savefig(allfigs_path(join(figure_folder, filename + ".pdf")))
    p.savefig(allfigs_path(join(figure_folder, filename + ".png")))
    # raw_input("Next plot?")
    p.close()

# Save table of parameters
co_param_df = DataFrame(co_fit_vals, index=g_CO.param_names)
co_param_df.to_latex(alltables_path("co_gaussian_totalprof_fits.tex"))
co_param_df.to_csv(iram_co21_14B088_data_path("tables/co_gaussian_totalprof_fits.csv",
                                              no_check=True))

# Per radial bin spectra
Nrows = 4
Ncols = 4

twocolumn_figure(font_scale=1.3)
p.figure(1, figsize=(8.4, 11)).clf()

fig, ax = p.subplots(Nrows, Ncols,
                     sharex=True,
                     sharey=True, num=1)

p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.text(0.5, 0.04, 'Velocity (km/s)', ha='center')
fig.text(0.04, 0.5, 'Normalized Intensity', va='center', rotation='vertical')


for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    r, c = np.unravel_index(ctr, (Nrows, Ncols))

    norm_hi = (total_spectrum_hi_radial[ctr] /
               total_spectrum_hi_radial[ctr].max()).value

    norm_co = (total_spectrum_co_radial[ctr] /
               total_spectrum_co_radial[ctr].max()).value

    ax[r, c].plot(hi_cube.spectral_axis.to(u.km / u.s).value,
                  norm_hi,
                  'b-', drawstyle='steps-mid', label="HI", alpha=0.7)
    # There's a 1 channel offset from my rotation subtraction in the cube
    ax[r, c].plot(co_cube.spectral_axis.to(u.km / u.s).value,
                  norm_co,
                  'g--', drawstyle='steps-mid', label="CO(2-1)", alpha=0.7)
    ax[r, c].set_ylim([-0.02, 1.1])
    ax[r, c].set_xlim([-120, 100])

    ax[r, c].annotate("{0} to {1}".format(r0.to(u.kpc).value, r1.to(u.kpc)),
                      xy=(-108, 0.65),
                      color='k', fontsize=8,
                      bbox={"boxstyle": "square", "facecolor": "w"})

    if ctr == 0:
        ax[r, c].legend(loc='upper left', frameon=True, prop={"size": 8})
    ax[r, c].grid()

for r in range(Nrows):
    for c in range(Ncols):
        if r == Nrows - 1:
            ax[r, c].set_xticklabels(ax[r, c].xaxis.get_majorticklabels(),
                                     rotation=45)

fig.savefig(allfigs_path(join(figure_folder, "total_profile_velocity_rotsub_hi_co_radial.pdf")))
fig.savefig(allfigs_path(join(figure_folder, "total_profile_velocity_rotsub_hi_co_radial.png")))

p.close()

# How do the model parameters change with radius?

g_CO_init = models.Gaussian1D(amplitude=1., mean=0., stddev=9.)
g_HI_init = models.Gaussian1D(amplitude=1., mean=0., stddev=10.)
# g_HI_init = models.Gaussian1D(amplitude=0.75, mean=0., stddev=5.) + \
#     models.Gaussian1D(amplitude=0.25, mean=0.0, stddev=13.)
# g_HI_init.mean_1.tied = tie_mean
# g_HI_init.amplitude_0.bounds = (0.5, 1.0)
# g_HI_init.amplitude_1.bounds = (0.0, 0.5)

hi_params = {}
co_params = {}

labels = ["rotsub", "centsub", "peaksub"]

for sub in labels:
    for name in g_HI_init.param_names:
        # Skip the tied mean
        if name == "mean_1":
            continue
        par_name = "{0}_{1}".format(sub, name)
        par_error = "{}_stderr".format(par_name)

        hi_params[par_name] = np.zeros_like(inneredge.value)
        hi_params[par_error] = np.zeros_like(inneredge.value)

    for name in g_CO_init.param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_error = "{}_stderr".format(par_name)

        co_params[par_name] = np.zeros_like(inneredge.value)
        co_params[par_error] = np.zeros_like(inneredge.value)


for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    hi_spectra = [total_spectrum_hi_radial[ctr],
                  total_spectrum_hi_radial_cent[ctr],
                  total_spectrum_hi_radial_peakvel[ctr]]

    for spectrum, label in zip(hi_spectra, labels):

        fit_g = fitting.LevMarLSQFitter()

        vels = hi_cube.spectral_axis.to(u.km / u.s).value
        norm_intens = (spectrum / spectrum.max()).value
        g_HI = fit_g(g_HI_init, vels, norm_intens, maxiter=1000)

        cov = fit_g.fit_info['param_cov']
        if cov is None:
            raise Exception("No covariance matrix")

        idx_corr = 0
        for idx, name in enumerate(g_HI.param_names):
            if name == "mean_1":
                idx_corr = 1
                continue
            par_name = "{0}_{1}".format(label, name)
            hi_params[par_name][ctr] = g_HI.parameters[idx]
            hi_params["{}_stderr".format(par_name)][ctr] = \
                np.sqrt(cov[idx - idx_corr, idx - idx_corr])

    co_spectra = [total_spectrum_co_radial[ctr],
                  total_spectrum_co_radial_cent[ctr],
                  total_spectrum_co_radial_peakvel[ctr]]

    for spectrum, label in zip(co_spectra, labels):

        fit_g = fitting.LevMarLSQFitter()

        co_vels = co_cube.spectral_axis.to(u.km / u.s).value
        norm_co_intens = (spectrum / spectrum.max()).value
        g_CO = fit_g_co(g_CO_init, co_vels, norm_co_intens, maxiter=1000)

        cov = fit_g_co.fit_info['param_cov']
        if cov is None:
            raise Exception("No covariance matrix")

        for idx, name in enumerate(g_CO.param_names):
            par_name = "{0}_{1}".format(label, name)
            co_params[par_name][ctr] = g_CO.parameters[idx]
            co_params["{}_stderr".format(par_name)][ctr] = \
                np.sqrt(cov[idx, idx])

bin_names = ["{}-{}".format(r0.value, r1)
             for r0, r1 in zip(inneredge, outeredge)]

co_radial_fits = DataFrame(co_params, index=bin_names)
hi_radial_fits = DataFrame(hi_params, index=bin_names)

co_radial_fits.to_latex(alltables_path("co_gaussian_totalprof_fits_radial.tex"))
co_radial_fits.to_csv(iram_co21_14B088_data_path("tables/co_gaussian_totalprof_fits_radial.csv",
                                                 no_check=True))

hi_radial_fits.to_latex(alltables_path("hi_gaussian_totalprof_fits_radial.tex"))
hi_radial_fits.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits_radial.csv",
                                             no_check=True))

# Plot comparisons of these fits
bin_cents = (outeredge - dr / 2.).to(u.kpc).value

hi_velres = \
  (hi_cube.spectral_axis[1] - hi_cube.spectral_axis[0]).to(u.km / u.s).value
co_velres = \
  (co_cube.spectral_axis[1] - co_cube.spectral_axis[0]).to(u.km / u.s).value

hi_width_error = lambda val: np.sqrt(val**2 + hi_velres**2)
co_width_error = lambda val: np.sqrt(val**2 + co_velres**2)

twocolumn_twopanel_figure()

fig, ax = p.subplots(1, 3, sharey=True)

ax[0].errorbar(bin_cents, hi_params["rotsub_stddev"],
               yerr=hi_width_error(hi_params["rotsub_stddev_stderr"]),
               color='b', label='HI',
               drawstyle='steps-mid')
ax[0].errorbar(bin_cents, co_params["rotsub_stddev"],
               yerr=co_width_error(co_params["rotsub_stddev_stderr"]),
               color='g', linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[0].legend(loc='lower left', frameon=True)
ax[0].grid()
ax[0].set_ylim([0.25, 15])
ax[0].text(1.3, 13.5, "Rotation subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel("Fitted Gaussian Width (km/s)")

ax[1].errorbar(bin_cents, hi_params["centsub_stddev"],
               yerr=hi_width_error(hi_params["centsub_stddev_stderr"]),
               color='b', label='HI',
               drawstyle='steps-mid')
ax[1].errorbar(bin_cents, co_params["centsub_stddev"],
               yerr=co_width_error(co_params["centsub_stddev_stderr"]),
               color='g', linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].text(1.3, 13.5, "Centroid subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[2].errorbar(bin_cents, hi_params["peaksub_stddev"],
               yerr=hi_width_error(hi_params["peaksub_stddev_stderr"]),
               color='b', label='HI',
               drawstyle='steps-mid')
ax[2].errorbar(bin_cents, co_params["peaksub_stddev"],
               yerr=co_width_error(co_params["peaksub_stddev_stderr"]),
               color='g', linestyle='--', label='CO(2-1)',
               drawstyle='steps-mid')
ax[2].grid()
ax[2].set_xlabel("Radius (kpc)")
ax[2].text(1.3, 13.5, "Peak Vel. subtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
p.tight_layout()
p.subplots_adjust(hspace=0.05,
                  wspace=0.05)

fig.savefig(allfigs_path(join(figure_folder, "total_profile_radial_widths_HI_CO21.pdf")))
fig.savefig(allfigs_path(join(figure_folder, "total_profile_radial_widths_HI_CO21.png")))

p.close()

default_figure()
