
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u


def radial_profile(gal, cube, dr=100 * u.pc, mom0=None,
                   max_rad=10 * u.kpc,
                   mass_conversion=0.019 * u.M_sun / (u.K * u.km / u.s),
                   restfreq=1.414 * u.GHz):
    if mom0 is None:
        mom0 = cube.moment0()
    radius = gal.radius(header=cube.header).to(u.kpc).value

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
        sdprof[ctr] = np.nansum(mom0[idx].value) / \
            np.sum(np.isfinite(mom0[idx].value))
        sdprof_sigma[ctr] = \
            np.sqrt(np.nansum((mom0[idx].value - sdprof[ctr])**2.) / \
                    np.sum(np.isfinite(mom0[idx].value)))
        radprof[ctr] = np.nanmean(radius[idx])

    # Re-apply some units
    radprof = radprof * u.kpc
    sdprof = sdprof * mom0.unit
    sdprof = sdprof_sigma * mom0.unit

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
    # Testing with a subcube to make things quicker
    # cube = cube[400:600]

    # Create a radial profile of HI
    rs, sd, sd_sigma = radial_profile(g, cube)

    # Add in creating a radial profile of the archival and the Arecibo
    # Also CO? I guess these should be on the same grid/resolution...

    p.errorbar(rs, sd, "bD-", yerr=sd_sigma)
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ K$^{-1}$km$^{-1}$s)")
    p.xlabel(r"R (kpc)")
