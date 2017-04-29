
from spectral_cube import SpectralCube
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.utils.console import ProgressBar


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
    from astropy.coordinates import Angle

    from paths import (fourteenB_HI_data_path, arecibo_HI_data_path,
                       c_hi_analysispath, paper1_figures_path,
                       data_path, paper1_tables_path)

    from constants import (lwidth_name, rotsub_cube_name,
                           rotsub_mask_name, hi_freq,
                           centroidsub_cube_name, centroidsub_mask_name,
                           peakvelsub_cube_name, peakvelssub_mask_name)

    from galaxy_params import gal
    from plotting_styles import default_figure, onecolumn_figure

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
    hi_pas = gal.position_angles(header=hi_cube.header)

    dr = 100 * u.pc
    max_radius = (8.0 * u.kpc).to(u.pc)

    nbins = np.int(np.floor(max_radius / dr))
    inneredge = np.linspace(0, max_radius - dr, nbins)
    outeredge = np.linspace(dr, max_radius, nbins)

    total_spectrum_hi_radial = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_n = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_s = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_cent = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_cent_n = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_cent_s = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_peakvel = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_peakvel_n = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
    total_spectrum_hi_radial_peakvel_s = \
        np.zeros((inneredge.size, hi_cube.shape[0])) * u.K

    # Make the N and S masks
    pa_bounds_n = Angle([0.5 * np.pi * u.rad, -0.5 * np.pi * u.rad])
    pa_bounds_n = pa_bounds_n.wrap_at(0.5 * np.pi * u.rad)
    hi_pa_mask_n = np.logical_and(hi_pas.wrap_at(0.5 * np.pi * u.rad) >= pa_bounds_n[0],
                                  hi_pas.wrap_at(0.5 * np.pi * u.rad) < pa_bounds_n[1])

    pa_bounds_s = Angle([-0.5 * np.pi * u.rad, 0.5 * np.pi * u.rad])
    hi_pa_mask_s = np.logical_and(hi_pas >= pa_bounds_s[0],
                                  hi_pas < pa_bounds_s[1])

    for ctr, (r0, r1) in enumerate(zip(inneredge,
                                       outeredge)):

        print("On bin {} to {}".format(r0.value, r1))

        hi_rad_mask = np.logical_and(hi_radius >= r0,
                                     hi_radius < r1)

        hi_mask_n = np.logical_and(hi_rad_mask, hi_pa_mask_n)
        hi_mask_s = np.logical_and(hi_rad_mask, hi_pa_mask_s)

        total_spectrum_hi_radial_n[ctr] = \
            total_profile(hi_cube, hi_mask_n).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))
        total_spectrum_hi_radial_s[ctr] = \
            total_profile(hi_cube, hi_mask_s).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

        total_spectrum_hi_radial[ctr] = \
            total_spectrum_hi_radial_n[ctr] + total_spectrum_hi_radial_s[ctr]

        total_spectrum_hi_radial_cent_n[ctr] = \
            total_profile(hi_cube_cent, hi_mask_n).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))
        total_spectrum_hi_radial_cent_s[ctr] = \
            total_profile(hi_cube_cent, hi_mask_s).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

        total_spectrum_hi_radial_cent[ctr] = \
            total_spectrum_hi_radial_cent_n[ctr] + \
            total_spectrum_hi_radial_cent_s[ctr]

        total_spectrum_hi_radial_peakvel_n[ctr] = \
            total_profile(hi_cube_peakvel, hi_mask_n).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))
        total_spectrum_hi_radial_peakvel_s[ctr] = \
            total_profile(hi_cube_peakvel, hi_mask_s).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

        total_spectrum_hi_radial_peakvel[ctr] = \
            total_spectrum_hi_radial_peakvel_n[ctr] + \
            total_spectrum_hi_radial_peakvel_s[ctr]


    # We'll make mock SpectralCubes from these so it's easy to calculate
    # moments
    rot_stack_n = SpectralCube(data=total_spectrum_hi_radial_n.T.reshape((1178, 80, 1)),
                               wcs=hi_cube.wcs)
    rot_stack_s = SpectralCube(data=total_spectrum_hi_radial_s.T.reshape((1178, 80, 1)),
                               wcs=hi_cube.wcs)
    rot_stack = SpectralCube(data=total_spectrum_hi_radial.T.reshape((1178, 80, 1)),
                             wcs=hi_cube.wcs)

    cent_stack_n = SpectralCube(data=total_spectrum_hi_radial_cent_n.T.reshape((1178, 80, 1)),
                                wcs=hi_cube.wcs)
    cent_stack_s = SpectralCube(data=total_spectrum_hi_radial_cent_s.T.reshape((1178, 80, 1)),
                                wcs=hi_cube.wcs)
    cent_stack = SpectralCube(data=total_spectrum_hi_radial_cent.T.reshape((1178, 80, 1)),
                              wcs=hi_cube.wcs)

    peakvel_stack_n = SpectralCube(data=total_spectrum_hi_radial_peakvel_n.T.reshape((1178, 80, 1)),
                                   wcs=hi_cube.wcs)
    peakvel_stack_s = SpectralCube(data=total_spectrum_hi_radial_peakvel_s.T.reshape((1178, 80, 1)),
                                   wcs=hi_cube.wcs)
    peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((1178, 80, 1)),
                                 wcs=hi_cube.wcs)

    # Now save all of these for future use.
    wstring = "{0}{1}".format(int(dr.value), dr.unit)
    rot_stack_n.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_north_{}.fits".format(wstring),
                                             no_check=True), overwrite=True)
    rot_stack_s.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_south_{}.fits".format(wstring),
                                             no_check=True), overwrite=True)
    rot_stack.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_{}.fits".format(wstring),
                                           no_check=True), overwrite=True)

    cent_stack_n.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_north_{}.fits".format(wstring),
                                              no_check=True), overwrite=True)
    cent_stack_s.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_south_{}.fits".format(wstring),
                                              no_check=True), overwrite=True)
    cent_stack.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_{}.fits".format(wstring),
                                            no_check=True), overwrite=True)

    peakvel_stack_n.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring),
                                                 no_check=True), overwrite=True)
    peakvel_stack_s.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring),
                                                 no_check=True), overwrite=True)
    peakvel_stack.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring),
                                               no_check=True), overwrite=True)

    rot_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_north_100pc.fits"))
    rot_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_south_100pc.fits"))
    rot_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_100pc.fits"))

    cent_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_north_100pc.fits"))
    cent_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_south_100pc.fits"))
    cent_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_100pc.fits"))

    peakvel_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_north_100pc.fits"))
    peakvel_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_south_100pc.fits"))
    peakvel_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_100pc.fits"))


    # Finally, fit Gaussian models and save the fit results

    g_HI_init = models.Gaussian1D(amplitude=1., mean=0., stddev=10.)

    hi_params = {}
    labels = ["rotsub", "rotsub_n", "rotsub_s", "centsub", "centsub_n",
              "centsub_s", "peaksub", "peaksub_n", "peaksub_s"]

    for sub in labels:
        for name in g_HI_init.param_names:
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
                      peakvel_stack_s[:, ctr, 0]]

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

    # Add stderr in quadrature with the channel width
    hi_velres = \
        (hi_cube.spectral_axis[1] -
         hi_cube.spectral_axis[0]).to(u.km / u.s).value

    # Add the velocity width of the channel in quadrature
    for col in hi_params.keys():
        if "mean_stderr" in col or "stddev_stderr" in col:
            hi_params[col + "_w_chanwidth"] = np.sqrt(hi_params[col]**2 +
                                                      hi_velres**2)

    hi_radial_fits = DataFrame(hi_params, index=bin_names)

    bin_string = "{0}{1}".format(int(dr.value), dr.unit)
    hi_radial_fits.to_latex(paper1_tables_path("hi_gaussian_totalprof_fits_radial_{}.tex".format(bin_string)))
    hi_radial_fits.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits_radial_{}.csv".format(bin_string),
                                                 no_check=True))
