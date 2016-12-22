
from astropy.utils.console import ProgressBar
import astropy.units as u
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import average_beams
import numpy as np
import matplotlib.pyplot as p
from pandas import DataFrame
from astropy.modeling import models, fitting
from astropy.io import fits

from analysis.paths import (fourteenB_HI_data_path, iram_co21_data_path,
                            paper1_figures_path, paper1_tables_path)
from analysis.constants import (rotsub_cube_name, rotsub_mask_name,
                                co21_mass_conversion, hi_freq,
                                centroidsub_cube_name, centroidsub_mask_name)
from analysis.galaxy_params import gal

'''
Create profiles of HI and CO after rotation subtraction.
'''


def total_profile(cube, spatial_mask=None, verbose=True):
    '''
    Create the total profile over a region in a given spatial mask.
    '''

    if verbose:
        chan_iter = ProgressBar(cube.shape[0])
    else:
        chan_iter = range(cube.shape[0])

    total_profile = np.zeros((cube.shape[0],)) * cube.unit

    for chan in chan_iter:
        if spatial_mask is not None:
            plane = cube[chan] * spatial_mask
        else:
            plane = cube[chan]
        total_profile[chan] = np.nansum(plane)

    return total_profile


co_cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.rotsub.fits"))
co_cube_cent = \
    SpectralCube.read(iram_co21_data_path("m33.co21_iram.hi_centroid_corrected.fits"))

hi_cube = SpectralCube.read(fourteenB_HI_data_path(rotsub_cube_name))
hi_mask = fits.open(fourteenB_HI_data_path(rotsub_mask_name))[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

hi_cube_cent = SpectralCube.read(fourteenB_HI_data_path(centroidsub_cube_name))
hi_mask_cent = fits.open(fourteenB_HI_data_path(centroidsub_mask_name))[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)


hi_beam = average_beams(hi_cube.beams)

hi_radius = gal.radius(header=hi_cube.header)
co_radius = gal.radius(header=co_cube.header)

# Perform the same analysis split up into radial bins
dr = 500 * u.pc

max_radius = (6.0 * u.kpc).to(u.pc)

nbins = np.int(np.floor(max_radius / dr))

inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)

total_spectrum_hi_radial = np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_co_radial = np.zeros((inneredge.size, co_cube.shape[0])) * u.K

total_spectrum_hi_radial_cent = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_co_radial_cent = \
    np.zeros((inneredge.size, co_cube.shape[0])) * u.K

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {}".format(r0.value, r1))

    hi_mask = np.logical_and(hi_radius >= r0,
                             hi_radius < r1)
    co_mask = np.logical_and(co_radius >= r0,
                             co_radius < r1)

    total_spectrum_hi_radial[ctr] = \
        total_profile(hi_cube, hi_mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_radial_cent[ctr] = \
        total_profile(hi_cube_cent, hi_mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_co_radial[ctr] = total_profile(co_cube, co_mask)

    total_spectrum_co_radial_cent[ctr] = total_profile(co_cube_cent, co_mask)

# Need to get portions of HI emission beyond 6 kpc.
total_spectrum_hi = \
    total_profile(hi_cube).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

total_spectrum_hi_cent = \
    total_profile(hi_cube_cent).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

# Significant CO emission is limited to within about 6 kpc
total_spectrum_co = total_spectrum_co_radial.sum(0)

total_spectrum_co_cent = total_spectrum_co_radial_cent.sum(0)


# Plot the profiles.
fig, ax = p.subplots(1, 2, sharey=True)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)
ax[0].plot(hi_cube.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi / total_spectrum_hi.max()).value,
           'b-', drawstyle='steps-mid', label="HI")
# There's a 1 channel offset from my rotation subtraction in the cube
ax[0].plot(co_cube.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co / total_spectrum_co.max()).value,
           'g--', drawstyle='steps-mid', label="CO(2-1)")
ax[0].set_xlabel("Velocity (km/s)")
ax[0].set_ylabel("Normalized Total Intensity")
ax[0].set_title("Rotation subtracted")
ax[0].set_ylim([-0.02, 1.1])
ax[0].set_xlim([-50, 50])
ax[0].grid()
ax[0].legend()

ax[1].plot(hi_cube_cent.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_hi_cent / total_spectrum_hi_cent.max()).value,
           'b-', drawstyle='steps-mid')
# There's a 1 channel offset from my rotation subtraction in the cube
ax[1].plot(co_cube.spectral_axis.to(u.km / u.s).value,
           (total_spectrum_co_cent / total_spectrum_co_cent.max()).value,
           'g--', drawstyle='steps-mid')
ax[1].set_title("Centroid subtracted")
ax[1].set_xlabel("Velocity (km/s)")
ax[1].set_ylim([-0.02, 1.1])
ax[1].set_xlim([-50, 50])
ax[1].grid()

p.draw()

fig.savefig(paper1_figures_path("total_profile_corrected_velocity_HI_CO21.pdf"))
fig.savefig(paper1_figures_path("total_profile_corrected_velocity_HI_CO21.png"))

p.clf()
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

vels = hi_cube.spectral_axis.to(u.km / u.s).value
norm_intens = (total_spectrum_hi / total_spectrum_hi.max()).value
g_HI = fit_g(g_HI_init, vels, norm_intens)

# The covariance matrix is hidden away... tricksy
cov = fit_g.fit_info['param_cov']
parnames = [n for n in g_HI.param_names if n not in ['mean_1']]
parvals = [v for (n, v) in zip(g_HI.param_names, g_HI.parameters)
           if n in parnames]

# Note that the statistical errors on the mean are too small.

p.plot(vels, norm_intens, 'b-', drawstyle='steps-mid')
p.plot(vels, g_HI(vels), 'k:', label="Total Fit")
p.plot(vels, g_HI["None_0"](vels), 'g--', label="Narrow Component")
p.plot(vels, g_HI["None_1"](vels), 'm-.', label="Wide Component")
p.xlabel("Velocity (km/s)")
p.ylabel("Total Normalized Intensity")
p.xlim([-100, 100])
p.legend()
p.ylim([-0.1, 1.1])
p.grid()
p.draw()

p.savefig(paper1_figures_path("total_profile_corrected_velocity_rotsub_hi_fit.pdf"))
p.savefig(paper1_figures_path("total_profile_corrected_velocity_rotsub_hi_fit.png"))
# raw_input("Next plot?")
p.clf()

fit_g = fitting.LevMarLSQFitter()

vels = hi_cube_cent.spectral_axis.to(u.km / u.s).value
norm_intens = (total_spectrum_hi_cent / total_spectrum_hi_cent.max()).value
g_HI = fit_g(g_HI_init, vels, norm_intens)

# The covariance matrix is hidden away... tricksy
cov_cent = fit_g.fit_info['param_cov']
parnames_cent = [n for n in g_HI.param_names if n not in ['mean_1']]
parvals_cent = [v for (n, v) in zip(g_HI.param_names, g_HI.parameters)
                if n in parnames]


# Save parameter table
hi_param_df = DataFrame({"Rot Sub Params": parvals,
                         "Rot Sub Errors": [np.sqrt(cov[i, i]) for i in
                                            range(cov.shape[0])],
                         "Cent Sub Params": parvals_cent,
                         "Cent Sub Errors": [np.sqrt(cov_cent[i, i]) for i in
                                             range(cov_cent.shape[0])]},
                        index=parnames)
hi_param_df.to_latex(paper1_tables_path("hi_gaussian_totalprof_fits.tex"))
hi_param_df.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits.csv",
                                          no_check=True))

# CO model

g_CO_init = models.Gaussian1D(amplitude=1., mean=0., stddev=9.)

fit_g_co = fitting.LevMarLSQFitter()

co_vels = co_cube.spectral_axis.to(u.km / u.s).value
norm_co_intens = total_spectrum_co / total_spectrum_co.max()
# norm_co_intens = np.roll((total_spectrum_co / total_spectrum_co.max()).value,
#                          -1)
g_CO = fit_g_co(g_CO_init, co_vels, norm_co_intens)

cov = fit_g_co.fit_info['param_cov']

# Better sampling for plotting
more_vels = np.arange(co_vels.min(), co_vels.max(), 0.5)

p.plot(co_vels, norm_co_intens, 'b-', drawstyle='steps-mid')
p.plot(more_vels, g_CO(more_vels), 'k--', label="Total Fit")
p.xlabel("Velocity (km/s)")
p.ylabel("Total Normalized Intensity")
p.xlim([-100, 100])
p.ylim([-0.1, 1.1])
p.grid()
p.draw()

p.savefig(paper1_figures_path("total_profile_corrected_velocity_co21_fit.pdf"))
p.savefig(paper1_figures_path("total_profile_corrected_velocity_co21_fit.png"))

p.close()


fit_g_co = fitting.LevMarLSQFitter()

co_vels = co_cube.spectral_axis.to(u.km / u.s).value
norm_co_intens = total_spectrum_co_cent / total_spectrum_co_cent.max()
g_CO_cent = fit_g_co(g_CO_init, co_vels, norm_co_intens)

cov_cent = fit_g_co.fit_info['param_cov']

# Save table of parameters
co_param_df = DataFrame({"Rot Sub Params": g_CO.parameters,
                         "Rot Sub Errors": [np.sqrt(cov[i, i]) for i in
                                            range(cov.shape[0])],
                         "Cent Sub Params": g_CO_cent.parameters,
                         "Cent Sub Errors": [np.sqrt(cov_cent[i, i]) for i in
                                             range(cov_cent.shape[0])]},
                        index=g_CO.param_names)
co_param_df.to_latex(paper1_tables_path("co_gaussian_totalprof_fits.tex"))
co_param_df.to_csv(iram_co21_data_path("tables/co_gaussian_totalprof_fits.csv",
                                       no_check=True))

# Per radial bin spectra
Nrows = 4
Ncols = 3

p.figure(1, figsize=(12, 20)).clf()

fig, ax = p.subplots(Nrows, Ncols,
                     sharex=True,
                     sharey=True, num=1)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)

fig.text(0.5, 0.02, 'Velocity (km/s)', ha='center')
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
    ax[r, c].set_xlim([-110, 100])

    ax[r, c].annotate("{0} to {1}".format(r0.to(u.kpc).value, r1.to(u.kpc)),
                      xy=(-98, 0.65),
                      color='k',
                      fontsize=13)

    if ctr == 0:
        ax[r, c].legend(loc='upper left', fontsize=14)
    ax[r, c].grid()

for r in range(Nrows):
    for c in range(Ncols):
        if r == Nrows - 1:
            ax[r, c].set_xticklabels(ax[r, c].xaxis.get_majorticklabels(),
                                     rotation=45)

fig.savefig(paper1_figures_path("total_profile_velocity_rotsub_hi_co_radial.pdf"))
fig.savefig(paper1_figures_path("total_profile_velocity_rotsub_hi_co_radial.png"))

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

for sub in ["rotsub", "centsub"]:
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

    # HI rotsub
    fit_g = fitting.LevMarLSQFitter()

    vels = hi_cube.spectral_axis.to(u.km / u.s).value
    norm_intens = (total_spectrum_hi_radial[ctr] /
                   total_spectrum_hi_radial[ctr].max()).value
    g_HI = fit_g(g_HI_init, vels, norm_intens, maxiter=1000)

    cov = fit_g.fit_info['param_cov']
    if cov is None:
        raise Exception("No covariance matrix")

    idx_corr = 0
    for idx, name in enumerate(g_HI.param_names):
        if name == "mean_1":
            idx_corr = 1
            continue
        par_name = "rotsub_{}".format(name)
        hi_params[par_name][ctr] = g_HI.parameters[idx]
        hi_params["{}_stderr".format(par_name)][ctr] = \
            np.sqrt(cov[idx - idx_corr, idx - idx_corr])

    # HI centsub
    fit_g = fitting.LevMarLSQFitter()

    vels = hi_cube_cent.spectral_axis.to(u.km / u.s).value
    norm_intens = (total_spectrum_hi_radial_cent[ctr] /
                   total_spectrum_hi_radial_cent[ctr].max()).value
    g_HI = fit_g(g_HI_init, vels, norm_intens, maxiter=1000)

    cov = fit_g.fit_info['param_cov']
    if cov is None:
        raise Exception("No covariance matrix")

    idx_corr = 0
    for idx, name in enumerate(g_HI.param_names):
        if name == "mean_1":
            idx_corr = 1
            continue
        par_name = "centsub_{}".format(name)
        hi_params[par_name][ctr] = g_HI.parameters[idx]
        hi_params["{}_stderr".format(par_name)][ctr] = \
            np.sqrt(cov[idx - idx_corr, idx - idx_corr])

    # CO rotsub
    fit_g_co = fitting.LevMarLSQFitter()

    co_vels = co_cube.spectral_axis.to(u.km / u.s).value
    norm_co_intens = (total_spectrum_co_radial[ctr] /
                      total_spectrum_co_radial[ctr].max()).value
    g_CO = fit_g_co(g_CO_init, co_vels, norm_co_intens, maxiter=1000)

    cov = fit_g_co.fit_info['param_cov']
    if cov is None:
        raise Exception("No covariance matrix")

    for idx, name in enumerate(g_CO.param_names):
        par_name = "rotsub_{}".format(name)
        co_params[par_name][ctr] = g_CO.parameters[idx]
        co_params["{}_stderr".format(par_name)][ctr] = np.sqrt(cov[idx, idx])

    # CO centsub
    fit_g_co = fitting.LevMarLSQFitter()

    co_vels = co_cube_cent.spectral_axis.to(u.km / u.s).value
    norm_co_intens = (total_spectrum_co_radial_cent[ctr] /
                      total_spectrum_co_radial_cent[ctr].max()).value
    g_CO = fit_g_co(g_CO_init, co_vels, norm_co_intens, maxiter=1000)

    cov = fit_g_co.fit_info['param_cov']
    if cov is None:
        raise Exception("No covariance matrix")

    for idx, name in enumerate(g_CO.param_names):
        par_name = "centsub_{}".format(name)
        co_params[par_name][ctr] = g_CO.parameters[idx]
        co_params["{}_stderr".format(par_name)][ctr] = np.sqrt(cov[idx, idx])

bin_names = ["{}-{}".format(r0.value, r1)
             for r0, r1 in zip(inneredge, outeredge)]

co_radial_fits = DataFrame(co_params, index=bin_names)
hi_radial_fits = DataFrame(hi_params, index=bin_names)

co_radial_fits.to_latex(paper1_tables_path("co_gaussian_totalprof_fits_radial.tex"))
co_radial_fits.to_csv(iram_co21_data_path("tables/co_gaussian_totalprof_fits_radial.csv",
                                          no_check=True))

hi_radial_fits.to_latex(paper1_tables_path("hi_gaussian_totalprof_fits_radial.tex"))
hi_radial_fits.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits_radial.csv",
                                             no_check=True))

# Plot comparisons of these fits
bin_cents = (outeredge - dr / 2.).to(u.kpc).value

fig, ax = p.subplots(1, 2, sharey=True)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)

ax[0].errorbar(bin_cents, hi_params["rotsub_stddev"],
               yerr=hi_params["rotsub_stddev_stderr"],
               color='b', label='HI')
ax[0].errorbar(bin_cents, co_params["rotsub_stddev"],
               yerr=co_params["rotsub_stddev_stderr"],
               color='g', linestyle='--', label='CO(2-1)')
ax[0].legend(loc='lower left')
ax[0].grid()
ax[0].set_title("Rotation subtracted")
ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel("Gaussian Width (km/s)")

ax[1].errorbar(bin_cents, hi_params["centsub_stddev"],
               yerr=hi_params["centsub_stddev_stderr"],
               color='b', label='HI')
ax[1].errorbar(bin_cents, co_params["centsub_stddev"],
               yerr=co_params["centsub_stddev_stderr"],
               color='g', linestyle='--', label='CO(2-1)')
ax[1].legend(loc='lower left')
ax[1].grid()
ax[1].set_title("Centroid subtracted")
ax[1].set_xlabel("Radius (kpc)")

fig.savefig(paper1_figures_path("total_profile_radial_widths_HI_CO21.pdf"))
fig.savefig(paper1_figures_path("total_profile_radial_widths_HI_CO21.png"))

fig.close()
