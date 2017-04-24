
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.utils.console import ProgressBar


def radial_profile(gal, moment, header=None, dr=100 * u.pc,
                   max_rad=10 * u.kpc, pa_bounds=None, beam=None):
    '''
    Create a radial profile of a moment 2 array, optionally with limits on
    the angles used.

    Parameters
    ----------
    moment : spectral_cube.Projection

    pa_bounds : list of 2 Angles, optional
        When specified, limits the angles used when calculating the profile.
        e.g. pa_bounds=Angle([0.0*u.rad, np.pi*u.rad]))
    '''

    if moment.ndim != 2:
        raise ValueError("mom0 must be 2 dimensional.")

    if header is None:
        header = moment.header

    if beam is None:
        # See if its in the projection
        if "beam" in moment.meta:
            beam = moment.meta["beam"]
        else:
            print("No beam attached to the Projection.")

    if beam is not None:
        beam_pix = beam.sr.to(u.deg**2) / (header["CDELT2"] * u.deg)**2

    radius = gal.radius(header=header).to(u.kpc).value
    if pa_bounds is not None:
        # Check if they are angles
        if len(pa_bounds) != 2:
            raise IndexError("pa_bounds must contain 2 angles.")
        if not isinstance(pa_bounds, Angle):
            raise TypeError("pa_bounds must be an Angle.")

        # Return the array of PAs in the galaxy frame
        pas = gal.position_angles(header=header)

        # If the start angle is greater than the end, we need to wrap about
        # the discontinuity
        if pa_bounds[0] > pa_bounds[1]:
            initial_start = pa_bounds[0].copy()
            pa_bounds = pa_bounds.wrap_at(initial_start)
            pas = pas.wrap_at(initial_start)

    if max_rad is not None:
        max_rad = max_rad.to(u.kpc).value
    else:
        max_rad = radius.max() - dr

    dr = dr.to(u.kpc).value

    nbins = np.int(np.floor(max_rad / dr))

    inneredge = np.linspace(0, max_rad, nbins)
    outeredge = np.linspace(0 + dr, max_rad + dr, nbins)
    sdprof = np.zeros(nbins)
    sdprof_sigma = np.zeros(nbins)
    radprof = np.zeros(nbins)

    if pa_bounds is not None:
        pa_idx = np.logical_and(pas >= pa_bounds[0], pas < pa_bounds[1])

    for ctr, (r0, r1) in enumerate(zip(inneredge,
                                       outeredge)):

        idx = np.logical_and(radius >= r0, radius < r1)
        if pa_bounds is not None:
            idx = np.logical_and(idx, pa_idx)

        radprof[ctr] = np.nanmean(radius[idx])

        sdprof[ctr] = np.nansum(moment[idx].value) / \
            np.sum(np.isfinite(radius[idx]))
        sdprof_sigma[ctr] = \
            np.sqrt(np.nansum((moment[idx].value - sdprof[ctr])**2.) /
                    np.sum(np.isfinite(radius[idx])))
        if beam is not None:
            sdprof_sigma[ctr] /= \
                np.sqrt(np.sum(np.isfinite(radius[idx])) / beam_pix)

    # Re-apply some units
    radprof = radprof * u.kpc
    sdprof = sdprof * moment.unit
    sdprof_sigma = sdprof_sigma * moment.unit

    return radprof, sdprof, sdprof_sigma


def total_profile(cube, spatial_mask=None, verbose=True):
    '''
    Create the total profile over a region in a given spatial mask.
    '''

    posns = np.where(spatial_mask)

    if verbose:
        spec_iter = ProgressBar(posns[0].size)
    else:
        spec_iter = xrange(posns[0].size)

    total_profile = np.zeros((cube.shape[0],)) * cube.unit

    for i in spec_iter:
        y, x = posns[0][i], posns[1][i]

        spec = cube[:, y, x]
        mask_spec = cube.mask.include(view=(slice(None), y, x))
        valid = np.logical_and(np.isfinite(spec), mask_spec)

        total_profile[valid] += spec[valid]

    return total_profile


if __name__ == "__main__":
    import matplotlib.pyplot as p
    from spectral_cube.cube_utils import average_beams
    from astropy.modeling import models, fitting
    from pandas import DataFrame

    from paths import (fourteenB_HI_data_path, arecibo_HI_data_path,
                       c_hi_analysispath, paper1_figures_path,
                       data_path, paper1_tables_path)

    from constants import (lwidth_name, rotsub_cube_name,
                           rotsub_mask_name, hi_freq,
                           centroidsub_cube_name, centroidsub_mask_name,
                           peakvelsub_cube_name, peakvelssub_mask_name)

    from galaxy_params import gal
    from plotting_styles import default_figure, onecolumn_figure


    lwidth_hdu = fits.open(fourteenB_HI_data_path(lwidth_name))[0]
    lwidth = Projection.from_hdu(lwidth_hdu).to(u.km / u.s)

    # Create a radial profile of the HI vel disp out to 8 kpc.
    # Beyond 8 kpc, noise is dominant. It may be some reflection of the
    # warp, but I don't trust it right now.
    rs, sd, sd_sigma = radial_profile(gal, lwidth, max_rad=8 * u.kpc)

    onecolumn_figure(font_scale=1.)
    p.errorbar(rs.value, sd.value,
               yerr=sd_sigma.value, fmt="-", color="b",
               drawstyle='steps-mid')
    p.xlabel("Radius (kpc)")
    p.ylabel("HI Velocity Dispersion (km/s)")
    p.grid()
    p.tight_layout()
    p.draw()

    p.savefig(paper1_figures_path("hi_veldisp_profile.png"))
    p.savefig(paper1_figures_path("hi_veldisp_profile.pdf"))

    # raw_input("Next plot?")
    p.clf()

    # Create the North and South portions.
    rs_n, sd_n, sd_sigma_n = \
        radial_profile(gal, lwidth, max_rad=8 * u.kpc,
                       pa_bounds=Angle([0.5 * np.pi * u.rad,
                                        -0.5 * np.pi * u.rad]))
    sd_n = sd_n.to(u.km / u.s)
    sd_sigma_n = sd_sigma_n.to(u.km / u.s)
    rs_s, sd_s, sd_sigma_s = \
        radial_profile(gal, lwidth, max_rad=8 * u.kpc,
                       pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                        0.5 * np.pi * u.rad]))
    sd_s = sd_s.to(u.km / u.s)
    sd_sigma_s = sd_sigma_s.to(u.km / u.s)

    # p.plot(rs.value, sd.value, "k-.",
    #        drawstyle='steps-mid', label="Total")
    p.plot(rs_n.value, sd_n.value, "r-.", label="North",
           drawstyle='steps-mid')
    p.plot(rs_s.value, sd_s.value, "g--", label="South",
           drawstyle='steps-mid')
    p.errorbar(rs.value, sd.value,
               yerr=sd_sigma.value, fmt="-", color="b",
               drawstyle='steps-mid', label='Total')
    p.xlabel("Radius (kpc)")
    p.ylabel("HI Velocity Dispersion (km/s)")
    p.grid()
    p.legend()
    p.tight_layout()
    p.draw()

    p.savefig(paper1_figures_path("hi_veldisp_profile_n_s.png"))
    p.savefig(paper1_figures_path("hi_veldisp_profile_n_s.pdf"))

    # raw_input("Next plot?")
    p.close()

    # Now do line stacking with the same bin size

    hi_cube = SpectralCube.read(fourteenB_HI_data_path(rotsub_cube_name))
    hi_mask = fits.open(fourteenB_HI_data_path(rotsub_mask_name))[0]
    hi_cube = hi_cube.with_mask(hi_mask.data > 0)

    hi_cube_cent = SpectralCube.read(fourteenB_HI_data_path(centroidsub_cube_name))
    hi_mask_cent = fits.open(fourteenB_HI_data_path(centroidsub_mask_name))[0]
    hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

    hi_cube_peakvel = SpectralCube.read(fourteenB_HI_data_path(peakvelsub_cube_name))
    hi_mask_peakvel = fits.open(fourteenB_HI_data_path(peakvelssub_mask_name))[0]
    hi_cube_peakvel = hi_cube_cent.with_mask(hi_mask_peakvel.data > 0)

    hi_beam = average_beams(hi_cube.beams)

    hi_radius = gal.radius(header=hi_cube.header)

    dr = 100 * u.pc
    max_radius = (8.0 * u.kpc).to(u.pc)

    nbins = np.int(np.floor(max_radius / dr))
    inneredge = np.linspace(0, max_radius - dr, nbins)
    outeredge = np.linspace(dr, max_radius, nbins)

    total_spectrum_hi_radial = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_cent = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_peakvel = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K

    for ctr, (r0, r1) in enumerate(zip(inneredge,
                                       outeredge)):

        print("On bin {} to {}".format(r0.value, r1))

        hi_rad_mask = np.logical_and(hi_radius >= r0,
                                     hi_radius < r1)

        total_spectrum_hi_radial[ctr] = \
            total_profile(hi_cube, hi_rad_mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

        total_spectrum_hi_radial_cent[ctr] = \
            total_profile(hi_cube_cent, hi_rad_mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

        total_spectrum_hi_radial_peakvel[ctr] = \
            total_profile(hi_cube_peakvel, hi_rad_mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    g_HI_init = models.Gaussian1D(amplitude=1., mean=0., stddev=10.)

    hi_params = {}
    labels = ["rotsub", "centsub", "peaksub"]

    for sub in labels:
        for name in g_HI_init.param_names:
            par_name = "{0}_{1}".format(sub, name)
            par_error = "{}_stderr".format(par_name)

            hi_params[par_name] = np.zeros_like(inneredge.value)
            hi_params[par_error] = np.zeros_like(inneredge.value)

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

    bin_names = ["{}-{}".format(r0.value, r1)
                 for r0, r1 in zip(inneredge, outeredge)]

    bin_center = (inneredge + dr / 2.).to(u.kpc)
    hi_params["bin_center"] = bin_center

    hi_radial_fits = DataFrame(hi_params, index=bin_names)

    bin_string = "{0}{1}".format(int(dr.value), dr.unit)
    hi_radial_fits.to_latex(paper1_tables_path("hi_gaussian_totalprof_fits_radial_{}.tex".format(bin_string)))
    hi_radial_fits.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits_radial_{}.csv".format(bin_string),
                                                 no_check=True))

    hi_velres = \
        (hi_cube.spectral_axis[1] -
         hi_cube.spectral_axis[0]).to(u.km / u.s).value

    # Add the velocity width of the channel in quadrature
    hi_width_error = lambda val: np.sqrt(val**2 + hi_velres**2)

    onecolumn_figure(font_scale=1.)
    p.errorbar(rs.value, sd.value,
               yerr=sd_sigma.value, fmt="-",
               drawstyle='steps-mid', label='Averaged Line Width')
    p.xlabel("Radius (kpc)")
    p.ylabel("HI Velocity Dispersion (km/s)")

    # Now the stacked fits
    p.errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['rotsub_stddev'],
               yerr=hi_width_error(hi_radial_fits['rotsub_stddev_stderr']),
               fmt='D', label='Rotation Stacked Fit', alpha=0.5)
    p.errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['centsub_stddev'],
               yerr=hi_width_error(hi_radial_fits['centsub_stddev_stderr']),
               fmt='o', label='Centroid Stacked Fit', alpha=0.5)
    p.errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_stddev'],
               yerr=hi_width_error(hi_radial_fits['peaksub_stddev_stderr']),
               fmt='^', label='Peak Stacked Fit', alpha=0.5)
    p.grid()
    p.legend(frameon=True)
    p.tight_layout()

    p.savefig(paper1_figures_path("hi_veldisp_w_stackedfits.png"))
    p.savefig(paper1_figures_path("hi_veldisp_w_stackedfits.pdf"))

    p.close()

    default_figure()
