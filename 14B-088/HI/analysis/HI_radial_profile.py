
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.coordinates import Angle

hi_mass_conversion = 0.019 * (u.M_sun / u.pc ** 2) / (u.K * u.km / u.s)


def surfdens_radial_profile(gal, header=None, cube=None,
                            dr=100 * u.pc, mom0=None,
                            max_rad=10 * u.kpc,
                            weight_type="area",
                            mass_conversion=hi_mass_conversion,
                            restfreq=1.414 * u.GHz,
                            pa_bounds=None, beam=None):
    '''
    Create a radial profile, optionally with limits on the angles used.

    Parameters
    ----------
    weight_type : "area" or "mass"
        Return either the area weighted profile (Sum Sigma * dA / Sum dA) or
        the mass weighted profile (Sum Sigma^2 dA / Sum Sigma dA). See
        Leroy et al. (2013) for a thorough description.
    pa_bounds : list of 2 Angles, optional
        When specified, limits the angles used when calculating the profile.
        e.g. pa_bounds=Angle([0.0*u.rad, np.pi*u.rad]))
    '''

    mom0 = mom0.squeeze()

    if weight_type not in ["area", "mass"]:
        raise ValueError("weight_type must be 'area' or 'mass'.")

    if mom0.ndim != 2:
        raise ValueError("mom0 must be 2 dimensional.")

    if mom0 is None:
        if cube is None:
            raise ValueError("Must give cube when not given mom0")
        mom0 = cube.moment0()

    if header is None and cube is None:
        raise TypeError("Either a header or cube must be given.")

    if header is None:
        if cube is not None:
            header = cube.header
        elif mom0 is not None:
            header = mom0.header

    if beam is None:
        if cube is not None:
            beam = cube.beam
        elif mom0 is not None:
            beam = mom0.meta["beam"]

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

        if weight_type == "area":
            sdprof[ctr] = np.nansum(mom0[idx].value) / \
                np.sum(np.isfinite(radius[idx]))
            sdprof_sigma[ctr] = \
                np.sqrt(np.nansum((mom0[idx].value - sdprof[ctr])**2.) /
                        np.sum(np.isfinite(radius[idx])))
        else:
            # Has to be 'mass' since this is checked at the beginning
            sdprof[ctr] = np.nansum(np.power(mom0[idx].value, 2)) / \
                np.nansum(mom0[idx].value)

            # Now the std has weights of Sigma * dA, which should be
            # normalized
            weights = mom0[idx].value / np.nansum(mom0[idx].value)

            # No denominator, since the weights were normalize to unity.
            sdprof_sigma[ctr] = np.sqrt(np.nansum(weights * (mom0[idx].value - sdprof[ctr])**2))

        # Rescale the sigma based on the number of independent samples
        if beam is not None:
            sdprof_sigma[ctr] /= \
                np.sqrt(np.sum(np.isfinite(radius[idx])) / beam_pix)
        radprof[ctr] = np.nanmean(radius[idx])

    # Re-apply some units
    radprof = radprof * u.kpc
    sdprof = sdprof * mom0.unit  # / u.pc ** 2
    sdprof_sigma = sdprof_sigma * mom0.unit  # / u.pc ** 2

    # Correct for the los inclinations
    sdprof *= np.cos(gal.inclination)
    sdprof_sigma *= np.cos(gal.inclination)

    # If in Jy/bm, convert to K.
    if cube is not None:
        unit = cube.unit

    # if mom0 is not None:
    #     unit = mom0.unit

    if unit.is_equivalent(u.Jy):
        # The beam units are sort of implied
        sdprof = sdprof.to(u.Jy * u.km / u.s)
        sdprof_sigma = sdprof_sigma.to(u.Jy * u.km / u.s)
        # Now convert to K
        if restfreq is not None:
            sdprof = sdprof * beam.jtok(restfreq)
            sdprof_sigma = sdprof_sigma * beam.jtok(restfreq)
        else:
            sdprof = sdprof * beam.jtok(mom0.header["RESTFREQ"])
            sdprof_sigma = sdprof_sigma * \
                beam.jtok(mom0.header["RESTFREQ"])

        sdprof = sdprof / u.Jy
        sdprof_sigma = sdprof_sigma / u.Jy

    if mass_conversion is not None:
        sdprof = sdprof * mass_conversion
        sdprof_sigma = sdprof_sigma * mass_conversion

    return radprof, sdprof, sdprof_sigma


if __name__ == "__main__":
    from galaxies import Galaxy
    import matplotlib.pyplot as p
    from astropy.table import Table
    from spectral_cube.lower_dimensional_structures import Projection
    from astropy.wcs import WCS
    from radio_beam import Beam
    from astropy.io import fits

    import os

    from paths import (fourteenB_HI_data_path, arecibo_HI_data_path,
                       c_hi_analysispath, paper1_figures_path,
                       data_path)

    from constants import hi_freq, moment0_name

    g = Galaxy("M33")

    mom0_hdu = fits.open(fourteenB_HI_data_path(moment0_name))[0]
    mom0 = Projection(mom0_hdu.data, wcs=WCS(mom0_hdu.header), unit=u.km / u.s)
    mom0.meta["beam"] = Beam.from_fits_header(mom0_hdu.header)

    # Bin size in pc
    dr = 100 * u.pc

    # Create a radial profile of HI
    rs, sd, sd_sigma = surfdens_radial_profile(g, mom0=mom0, dr=dr,
                                               restfreq=hi_freq)
    rs_n, sd_n, sd_sigma_n = \
        surfdens_radial_profile(g, mom0=mom0,
                                pa_bounds=Angle([0.5 * np.pi * u.rad,
                                                -0.5 * np.pi * u.rad]),
                                dr=dr, restfreq=hi_freq)
    rs_s, sd_s, sd_sigma_s = \
        surfdens_radial_profile(g, mom0=mom0,
                                pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                                 0.5 * np.pi * u.rad]),
                                dr=dr, restfreq=hi_freq)

    # There is a ~1.5 global scaling factor b/w Arecibo and VLA + Arecibo
    # Not needed with the corrected cube!
    # scale_factor = 1.451
    scale_factor = 1.0

    # Add in creating a radial profile of the archival and the Arecibo
    # Also CO? I guess these should be on the same grid/resolution...

    # Arecibo
    # arecibo_file = arecibo_HI_data_path("14B-088_items/M33_14B-088_HI_model.fits")
    arecibo_file = arecibo_HI_data_path("M33only_jy_stokes_vrad.fits")
    arecibo_cube = \
        SpectralCube.read(arecibo_file)
    # Cut down to extents of the other cube
    # arecibo_cube = \
    #    arecibo_cube.subcube(xlo=cube.longitude_extrema[0],
    #                         xhi=cube.longitude_extrema[1],
    #                         yhi=cube.latitude_extrema[0],
    #                         ylo=cube.latitude_extrema[1])
    arecibo_mom0 = arecibo_cube.moment0()
    rs_arec, sd_arec, sd_sigma_arec = \
        surfdens_radial_profile(g, cube=arecibo_cube, mom0=arecibo_mom0, dr=dr,
                                restfreq=hi_freq)

    # Archival HI
    arch_vla_file = os.path.join(data_path, "VLA/AT0206/old_imaging/m33_hi.masked.fits")
    # arch_vla_file = os.path.join(data_path, "VLA/AT0206/imaging/M33_206_b_c_HI.fits")
    arch_cube = SpectralCube.read(arch_vla_file)

    arch_mom0 = arch_cube.moment0()
    rs_arch, sd_arch, sd_sigma_arch = \
        surfdens_radial_profile(g, cube=arch_cube, mom0=arch_mom0, dr=dr,
                                restfreq=hi_freq)

    p.ioff()

    # Show the total radial profile VLA and Arecibo
    p.errorbar(rs.value, sd.value / scale_factor,
               yerr=sd_sigma.value / scale_factor, fmt="-", color="b",
               label="VLA + Arecibo", drawstyle='steps-mid')
    p.plot(rs_arec.value, sd_arec.value, "g--", drawstyle='steps-mid',
           label="Arecibo")
    # p.errorbar(rs_arec.value, sd_arec.value, yerr=sd_sigma_arec.value,
    #            fmt="o--", color="g", label="Arecibo", drawstyle='steps-mid')
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    p.xlabel(r"R (kpc)")
    p.legend(loc='best')
    p.grid("on")
    p.savefig(paper1_figures_path("M33_Sigma_profile_w_Arecibo.pdf"))
    p.savefig(paper1_figures_path("M33_Sigma_profile_w_Arecibo.png"))
    p.clf()

    # W/ archival VLA
    p.errorbar(rs.value, sd.value / scale_factor,
               yerr=sd_sigma.value / scale_factor, fmt="-", color="b",
               label="VLA + Arecibo", drawstyle='steps-mid')
    p.plot(rs_arec.value, sd_arec.value, "g--", drawstyle='steps-mid',
           label="Arecibo")
    p.errorbar(rs_arch.value, sd_arch.value, yerr=sd_sigma_arec.value,
               color="r", fmt="-.", drawstyle='steps-mid',
               label="Archival VLA")
    # p.errorbar(rs_arec.value, sd_arec.value, yerr=sd_sigma_arec.value,
    #            fmt="o--", color="g", label="Arecibo", drawstyle='steps-mid')
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    p.xlabel(r"R (kpc)")
    p.legend(loc='best')
    p.grid("on")
    p.savefig(paper1_figures_path("M33_Sigma_profile_w_Arecibo_archival.pdf"))
    p.savefig(paper1_figures_path("M33_Sigma_profile_w_Arecibo_archival.png"))
    p.clf()

    # Show the north vs south profiles
    p.plot(rs.value, sd.value / scale_factor, "k-.",
           drawstyle='steps-mid', label="Total")
    p.plot(rs_n.value, sd_n.value / scale_factor, "b-", label="North",
           drawstyle='steps-mid')
    p.plot(rs_s.value, sd_s.value / scale_factor, "g--", label="South",
           drawstyle='steps-mid')
    # p.errorbar(rs_n.value, sd_n.value, yerr=sd_sigma_n.value, fmt="D-",
    #            color="b", label="North")
    # p.errorbar(rs_s.value, sd_s.value, yerr=sd_sigma_s.value, fmt="o-",
    #            color="g", label="South")
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    p.xlabel(r"R (kpc)")
    p.legend(loc='best')
    p.grid("on")
    p.savefig(paper1_figures_path("M33_Sigma_profile_N_S.pdf"))
    p.savefig(paper1_figures_path("M33_Sigma_profile_N_S.png"))
    p.clf()

    # Compare to the surface density profile in Corbelli
    corbelli_filename = \
        c_hi_analysispath("rotation_curves/corbelli_rotation_curve.csv")
    corbelli = Table.read(corbelli_filename)

    p.plot(rs.value, sd.value / scale_factor,
           linestyle="-", color="b",
           label="This work", drawstyle='steps-mid')
    p.plot(corbelli["R"][corbelli["R"] <= 10.0],
           corbelli["SigmaHI"][corbelli["R"] <= 10.0], "r--", drawstyle='steps-mid',
           label="Corbelli et al. (2014)")
    p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
    p.xlabel(r"R (kpc)")
    p.legend(loc='best')
    p.grid()
    p.savefig(paper1_figures_path("M33_Sigma_profile_w_Corbelli.pdf"))
    p.savefig(paper1_figures_path("M33_Sigma_profile_w_Corbelli.png"))
    p.clf()

    p.ion()
