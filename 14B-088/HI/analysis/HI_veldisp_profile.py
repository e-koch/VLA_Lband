
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


if __name__ == "__main__":
    import matplotlib.pyplot as p
    from spectral_cube.cube_utils import average_beams
    from astropy.modeling import models, fitting
    from pandas import DataFrame, read_csv

    from paths import (fourteenB_HI_data_path, arecibo_HI_data_path,
                       c_hi_analysispath, paper1_figures_path,
                       data_path, paper1_tables_path)

    from constants import (lwidth_name, rotsub_cube_name,
                           rotsub_mask_name, hi_freq,
                           centroidsub_cube_name, centroidsub_mask_name,
                           peakvelsub_cube_name, peakvelssub_mask_name)

    from galaxy_params import gal
    from plotting_styles import (default_figure, onecolumn_figure,
                                 twocolumn_twopanel_figure)


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

    # Now load in the line stacking fits with the same bin size

    hi_radial_fits = \
        read_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits_radial_100pc.csv"))

    onecolumn_figure(font_scale=1.)
    p.errorbar(rs.value, sd.value,
               yerr=sd_sigma.value, fmt="-",
               drawstyle='steps-mid', label='Averaged Line Width')
    p.xlabel("Radius (kpc)")
    p.ylabel("HI Velocity Dispersion (km/s)")

    # Now the stacked fits
    p.errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['rotsub_stddev'],
               yerr=hi_radial_fits['rotsub_stddev_stderr_w_chanwidth'],
               fmt='D', label='Rotation Stacked Fit', alpha=0.5)
    p.errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['centsub_stddev'],
               yerr=hi_radial_fits['centsub_stddev_stderr_w_chanwidth'],
               fmt='o', label='Centroid Stacked Fit', alpha=0.5)
    p.errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_stddev'],
               yerr=hi_radial_fits['peaksub_stddev_stderr_w_chanwidth'],
               fmt='^', label='Peak Stacked Fit', alpha=0.5)
    p.grid()
    p.legend(frameon=True)
    p.tight_layout()

    p.savefig(paper1_figures_path("hi_veldisp_w_stackedfits.png"))
    p.savefig(paper1_figures_path("hi_veldisp_w_stackedfits.pdf"))

    p.close()

    # Let's compare the line width from the second moment to the Gaussian width
    rot_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_100pc.fits"))
    cent_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_100pc.fits"))
    peakvel_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_100pc.fits"))

    twocolumn_twopanel_figure(font_scale=1.2)

    fig, ax = p.subplots(1, 3, sharey=True)
    ax[0].errorbar(hi_radial_fits['bin_center'],
                   hi_radial_fits['rotsub_stddev'],
                   yerr=hi_radial_fits['rotsub_stddev_stderr_w_chanwidth'],
                   fmt='D', label='Gaussian Fit', alpha=0.5)
    ax[0].plot(hi_radial_fits['bin_center'],
               rot_stack.linewidth_sigma().value / 1000.,
               label="Moment")
    ax[0].text(5, 11.5, "Rotation\nsubtracted",
               bbox={"boxstyle": "square", "facecolor": "w"})
    ax[0].legend(frameon=True, loc='lower right')
    ax[0].grid()
    # ax[0].set_xticklabels([])
    ax[0].set_ylabel("HI Velocity Dispersion (km/s)")
    ax[0].set_xlabel("Radius (kpc)")

    ax[1].errorbar(hi_radial_fits['bin_center'],
                   hi_radial_fits['centsub_stddev'],
                   yerr=hi_radial_fits['centsub_stddev_stderr_w_chanwidth'],
                   fmt='D', label='Gaussian Fit', alpha=0.5)
    ax[1].plot(hi_radial_fits['bin_center'],
               cent_stack.linewidth_sigma().value / 1000.,
               label="Moment")
    ax[1].text(5, 11.5, "Centroid\nsubtracted",
               bbox={"boxstyle": "square", "facecolor": "w"})
    ax[1].grid()
    ax[1].set_xlabel("Radius (kpc)")

    ax[2].errorbar(hi_radial_fits['bin_center'],
                   hi_radial_fits['peaksub_stddev'],
                   yerr=hi_radial_fits['peaksub_stddev_stderr_w_chanwidth'],
                   fmt='D', label='Gaussian Fit', alpha=0.5)
    ax[2].plot(hi_radial_fits['bin_center'],
               cent_stack.linewidth_sigma().value / 1000.,
               label="Moment")
    ax[2].text(5, 11.5, "Peak Vel.\nsubtracted",
               bbox={"boxstyle": "square", "facecolor": "w"})
    ax[2].grid()
    ax[2].set_xlabel("Radius (kpc)")

    p.tight_layout()

    fig.savefig(paper1_figures_path("hi_veldisp_avg_vs_stackedfits.png"))
    fig.savefig(paper1_figures_path("hi_veldisp_avg_vs_stackedfits.pdf"))
    p.close()

    default_figure()
