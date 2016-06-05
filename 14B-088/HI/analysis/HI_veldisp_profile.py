
from spectral_cube import SpectralCube
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle


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
    from galaxies import Galaxy
    import matplotlib.pyplot as p
    import os
    from astropy.io.fits import getdata

    save_lwidth = False

    direc = "/home/eric/MyRAID/M33/14B-088/HI/full_imaging/"

    cube_file = os.path.join(direc,
                             "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits")

    cube = \
        SpectralCube.read(cube_file)

    # Now load in the source mask
    cube = cube.with_mask(getdata(os.path.join(direc,
        "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub_source_mask.fits")).astype(bool))

    if save_lwidth:
        lwidth = cube.linewidth_sigma()
        lwidth.write(os.path.join(direc,
            "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.lwidth.fits"))
    else:
        load = fits.open(os.path.join(direc,
            "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.lwidth.fits")[0]
        lwidth = Projection(load.data, wcs=WCS(load.header))
        lwidth.meta["beam"] = cube.beam

    g = Galaxy("M33")

    # Create a radial profile of the HI vel disp out to 8 kpc.
    # Beyond 8 kpc, noise is dominant. It may be some reflection of the
    # warp, but I don't trust it right now.
    rs, sd, sd_sigma = radial_profile(g, lwidth, max_rad=8 * u.kpc)

    sd = sd.to(u.km / u.s)
    sd_sigma = sd_sigma.to(u.km / u.s)

    p.errorbar(rs.value, sd.value,
               yerr=sd_sigma.value, fmt="-", color="b",
               drawstyle='steps-mid')
    p.xlabel("R (kpc)")
    p.ylabel("HI Velocity Dispersion (km/s)")
    p.grid()
    p.draw()

    raw_input("Next plot?")
    p.clf()

    # Create the North and South portions.
    rs_n, sd_n, sd_sigma_n = \
        radial_profile(g, lwidth, max_rad=8 * u.kpc,
                       pa_bounds=Angle([0.5 * np.pi * u.rad,
                                        -0.5 * np.pi * u.rad]))
    sd_n = sd_n.to(u.km / u.s)
    sd_sigma_n = sd_sigma_n.to(u.km / u.s)
    rs_s, sd_s, sd_sigma_s = \
        radial_profile(g, lwidth, max_rad=8 * u.kpc,
                       pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                        0.5 * np.pi * u.rad]))
    sd_s = sd_s.to(u.km / u.s)
    sd_sigma_s = sd_sigma_s.to(u.km / u.s)

    p.plot(rs.value, sd.value, "k-.",
           drawstyle='steps-mid', label="Total")
    p.plot(rs_n.value, sd_n.value, "b-", label="North",
           drawstyle='steps-mid')
    p.plot(rs_s.value, sd_s.value, "g--", label="South",
           drawstyle='steps-mid')
    p.xlabel("R (kpc)")
    p.ylabel("HI Velocity Dispersion (km/s)")
    p.grid()
    p.legend()
    p.draw()
