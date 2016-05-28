
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.coordinates import Angle

hi_mass_conversion = 0.019 * (u.M_sun / u.pc **2) / (u.K * u.km / u.s)


def radial_profile(gal, cube, dr=100 * u.pc, mom0=None,
                   max_rad=10 * u.kpc,
                   mass_conversion=hi_mass_conversion,
                   restfreq=1.414 * u.GHz,
                   pa_bounds=None):
    '''
    Create a radial profile, optionally with limits on the angles used.

    Parameters
    ----------
    pa_bounds : list of 2 Angles, optional
        When specified, limits the angles used when calculating the profile.
        e.g. pa_bounds=Angle([0.0*u.rad, np.pi*u.rad]))
    '''

    if mom0 is None:
        mom0 = cube.moment0()
    radius = gal.radius(header=cube.header).to(u.kpc).value
    if pa_bounds is not None:
        # Check if they are angles
        if len(pa_bounds) != 2:
            raise IndexError("pa_bounds must contain 2 angles.")
        if not isinstance(pa_bounds, Angle):
            raise TypeError("pa_bounds must be an Angle.")

        # Return the array of PAs in the galaxy frame
        pas = gal.position_angles(header=cube.header)

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
        sdprof[ctr] = np.nansum(mom0[idx].value) / \
            np.sum(np.isfinite(radius[idx]))
        sdprof_sigma[ctr] = \
            np.sqrt(np.nansum((mom0[idx].value - sdprof[ctr])**2.) /
                    np.sum(np.isfinite(radius[idx])))
        radprof[ctr] = np.nanmean(radius[idx])

    # Re-apply some units
    radprof = radprof * u.kpc
    sdprof = sdprof * mom0.unit  # / u.pc ** 2
    sdprof_sigma = sdprof_sigma * mom0.unit  # / u.pc ** 2

    # Correct for the los inclinations
    sdprof *= np.cos(gal.inclination)
    sdprof_sigma *= np.cos(gal.inclination)

    # If in Jy/bm, convert to K.
    if cube.unit.is_equivalent(u.Jy):
        # The beam units are sort of implied
        sdprof = sdprof.to(u.Jy * u.km / u.s)
        sdprof_sigma = sdprof_sigma.to(u.Jy * u.km / u.s)
        # Now convert to K
        if restfreq is not None:
            sdprof = sdprof * mom0.meta["beam"].jtok(restfreq)
            sdprof_sigma = sdprof_sigma * mom0.meta["beam"].jtok(restfreq)
        else:
            sdprof = sdprof * mom0.meta["beam"].jtok(mom0.header["RESTFREQ"])
            sdprof_sigma = sdprof_sigma * \
                mom0.meta["beam"].jtok(mom0.header["RESTFREQ"])

        sdprof = sdprof / u.Jy
        sdprof_sigma = sdprof_sigma / u.Jy

    if mass_conversion is not None:
        sdprof = sdprof * mass_conversion
        sdprof_sigma = sdprof_sigma * mass_conversion

    return radprof, sdprof, sdprof_sigma

if __name__ == "__main__":
    from galaxies import Galaxy
    import matplotlib.pyplot as p
    import os

    direc = "/home/eric/MyRAID/M33/14B-088/HI/full_imaging/"
    # direc = "/home/eric/MyRAID/M33/14B-088/HI/combined_HI/"

    cube_file = os.path.join(direc,
                             "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits")
    # cube_file = os.path.join(direc,
    #                          "M33_14B-088_AT0206_HI.clean.image.pbcov_gt_0.3_masked.fits")

    cube = \
        SpectralCube.read(cube_file)

    g = Galaxy("M33")

    cube = cube.with_mask(cube > 1.8e-3 * u.Jy)

    mom0 = cube.moment0()

    # Create a radial profile of HI
    rs, sd, sd_sigma = radial_profile(g, cube, mom0=mom0)
    rs_n, sd_n, sd_sigma_n = \
        radial_profile(g, cube, mom0=mom0,
                       pa_bounds=Angle([0.5 * np.pi * u.rad,
                                        -0.5 * np.pi * u.rad]))
    rs_s, sd_s, sd_sigma_s = \
        radial_profile(g, cube, mom0=mom0,
                       pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                        0.5 * np.pi * u.rad]))

    # There is a ~1.5 global scaling factor b/w Arecibo and VLA + Arecibo
    scale_factor = 1.451

    # Add in creating a radial profile of the archival and the Arecibo
    # Also CO? I guess these should be on the same grid/resolution...

    # Arecibo
    arecibo_path = "/media/eric/Data_3//M33/Arecibo/14B-088_items/"
    arecibo_cube = \
        SpectralCube.read(os.path.join(arecibo_path,
                                       "M33_14B-088_HI_model.fits"))
    # Cut down to extents of the other cube
    arecibo_cube = \
        arecibo_cube.subcube(xlo=cube.longitude_extrema[0],
                             xhi=cube.longitude_extrema[1],
                             yhi=cube.latitude_extrema[0],
                             ylo=cube.latitude_extrema[1])
    arecibo_mom0 = arecibo_cube.moment0()
    rs_arec, sd_arec, sd_sigma_arec = \
        radial_profile(g, arecibo_cube, mom0=arecibo_mom0)


    p.ioff()

    # Show the total radial profile VLA and Arecibo
    p.errorbar(rs.value, sd.value / scale_factor,
               yerr=sd_sigma.value / scale_factor, fmt="D-", color="b",
               label="VLA + Arecibo")
    p.errorbar(rs_arec.value, sd_arec.value, yerr=sd_sigma_arec.value,
               fmt="o--", color="g", label="Arecibo")
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    p.xlabel(r"R (kpc)")
    p.legend(loc='best')
    p.grid("on")
    p.show()

    # Show the north vs south profiles
    p.plot(rs.value, sd.value / scale_factor, "k--")
    p.plot(rs_n.value, sd_n.value / scale_factor, "bD-", label="North")
    p.plot(rs_s.value, sd_s.value / scale_factor, "go-", label="South")
    # p.errorbar(rs_n.value, sd_n.value, yerr=sd_sigma_n.value, fmt="D-", color="b",
    #            label="North")
    # p.errorbar(rs_s.value, sd_s.value, yerr=sd_sigma_s.value, fmt="o-", color="g",
    #            label="South")
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    p.xlabel(r"R (kpc)")
    p.legend(loc='best')
    p.grid("on")
    p.show()

    p.ion()
