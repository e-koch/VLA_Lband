
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.coordinates import Angle

hi_mass_conversion = 0.019 * (u.M_sun / u.pc**2) / (u.K * u.km / u.s)


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
    for ctr, (r0, r1) in enumerate(zip(inneredge,
                                       outeredge)):

        idx = np.logical_and(radius >= r0, radius < r1)
        if pa_bounds is not None:
            pa_idx = np.logical_and(pas >= pa_bounds[0], pas < pa_bounds[1])
            idx = np.logical_and(idx, pa_idx)
        sdprof[ctr] = np.nansum(mom0[idx].value) / \
            np.sum(np.isfinite(mom0[idx].value))
        sdprof_sigma[ctr] = \
            np.sqrt(np.nansum((mom0[idx].value - sdprof[ctr])**2.) /
                    np.sum(np.isfinite(mom0[idx].value)))
        radprof[ctr] = np.nanmean(radius[idx])

    # Re-apply some units
    radprof = radprof * u.kpc
    sdprof = sdprof * mom0.unit
    sdprof_sigma = sdprof_sigma * mom0.unit

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

    cube_file = os.path.join(direc,
                             "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits")

    cube = \
        SpectralCube.read(cube_file)

    g = Galaxy("M33")

    cube = cube.with_mask(cube > 1.8e-3 * u.Jy)

    mom0 = cube.moment0()

    # Create a radial profile of HI
    rs, sd, sd_sigma = radial_profile(g, cube, mom0=mom0)

    # Add in creating a radial profile of the archival and the Arecibo
    # Also CO? I guess these should be on the same grid/resolution...

    p.errorbar(rs.value, sd.value, yerr=sd_sigma.value, fmt="D-", color="b")
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    p.xlabel(r"R (kpc)")
