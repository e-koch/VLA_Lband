
import numpy as np
from itertools import cycle
import scipy.optimize as opt

'''
Krumholz model for the fraction of H2 as a function of surface density and
metallicity.

Following Schruba+11, the output is in terms of the H2-HI fraction rather than
the the molecular fraction.
'''


def krumholz_ratio_model(Sigma, Z=0.5, c=1):
    '''

    Parameters
    ----------
    sigma : float or np.ndarray
        Surface Density in Msol/pc^2
    Z : float
        Metallicity in units of solar metallicity.
    c : float
        Clumping fraction. Expected to be near unity with a resolution of
        100 pc. c>=1.
    '''

    Sigma_comp = c * Sigma

    chi = 0.77 * (1 + 3.1 * np.power(Z, 0.365))

    s = np.log(1 + 0.6 * chi) / (0.04 * Sigma_comp * Z)

    delta = 0.0712 * np.power((0.1 / s) + 0.675, -2.8)

    frac = 1 - np.power(1 + np.power(0.75 * (s / (1 + delta)), -5.), -0.2)

    return frac / (1 - frac)


def alternate_krumholz_ratio_model(Sigma, psi=1.0, c=1, Z=0.1):
    '''

    Parameters
    ----------
    sigma : float or np.ndarray
        Surface Density in Msol/pc^2
    psi : float
        Dust-adjusted radiation field (unitless). Related to metallicity
        (among other things). At Z=1, psi=1.6, and at Z=0.1, psi=1.0
    c : float
        Clumping fraction. Expected to be near unity with a resolution of
        100 pc. c>=1.
    '''

    Sigma_comp = c * Sigma

    s = Sigma_comp * Z / float(psi)

    term1 = (s / 11.) ** 3 * ((125 + s) / (96 + s)) ** 3

    return np.power(1 + term1, 1 / 3.) - 1


def optimize_clump_factors(Sigma, R, Z=0.5):
    '''
    Solve for the clump factor needed to intersect each point.

    Parameters
    ----------
    Sigma : Quantity
        Surface densities.
    Z : float or array
        Metallicity values.
    '''

    if not isinstance(Z, np.ndarray):
        Z = cycle([Z])

    clump_values = []

    for sig, r, z in zip(Sigma, R, Z):

        def model(sigma, c):
            return krumholz_ratio_model(sigma, Z=z, c=c)

        popt, pcov = opt.curve_fit(model, sig, r)

        clump_values.append(popt[0])

    return np.array(clump_values)


if __name__ == "__main__":

    import matplotlib.pyplot as p
    from astropy import units as u
    from spectral_cube import SpectralCube
    import astropy.io.fits as fits
    from astropy.table import Table
    from reproject import reproject_interp
    from radio_beam import Beam
    from scipy.stats import binned_statistic

    from paths import (iram_co21_data_path, fourteenB_HI_data_path,
                       paper1_figures_path)
    from constants import (moment0_name,
                           co21_mass_conversion, hi_mass_conversion,
                           hi_freq)
    from galaxy_params import gal
    from plotting_styles import onecolumn_figure, default_figure

    # Bin radius should match whatever was used in co_radial_profile.py
    dr = 100.0 * u.pc

    # Read in the radial profiles table
    tab_name = "tables/co21_hi_radialprofiles_{}pc.fits".format(int(dr.value))
    try:
        tab = Table.read(fourteenB_HI_data_path(tab_name))
    except OSError:
        raise OSError("Table does not exist in the 14B-088 data path. "
                      "Run co_radial_profile.py first:"
                      " {}".format(tab_name))

    rs = u.Quantity(tab["Radius"])
    sd = u.Quantity(tab["CO_Sigma"])
    sd_sigma = u.Quantity(tab["CO_Sigma_std"])
    sd_hi = u.Quantity(tab["HI_Sigma"])
    sd_sigma_hi = u.Quantity(tab["HI_Sigma_std"])

    # Now plot their ratio against the total gas surface density
    gas_ratio = sd.value / sd_hi.value
    gas_ratio_sigma = \
        (gas_ratio *
         np.sqrt((sd_sigma / sd)**2 + (sd_sigma_hi / sd_hi)**2)).value
    log_gas_ratio_sigma = gas_ratio_sigma / (gas_ratio * np.log(10))
    total_sd = sd.value + sd_hi.value
    total_sd_sigma = \
        (total_sd *
         np.sqrt((sd_sigma / sd)**2 + (sd_sigma_hi / sd_hi)**2)).value

    # Overplot the Krumholz model with a few different clumping factors.
    # Theoretically, c -> 1 at a resolution of 100 pc. but I'm finding a better
    # match when c=4-6. The model is supposed to take metallicity into account,
    # but maybe the gradient is causing some issues? Schruba+11 finds c~2 for
    # their entire sample, with a lot of scatter
    sds = np.arange(1, 40, 0.2)

    onecolumn_figure(font_scale=1.0, fig_ratio=1.0)
    # p.semilogy(total_sd, gas_ratio, 'bD')
    p.errorbar(total_sd, np.log10(gas_ratio), yerr=log_gas_ratio_sigma,
               xerr=total_sd_sigma, color='b', alpha=0.6, fmt='D')
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=1, Z=0.5)), "r--",
           label="c=1, Z=0.5")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=4, Z=0.5)), "g-.",
           label="c=4, Z=0.5")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=4, Z=0.25)), "m.",
           label="c=4, Z=0.25")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=4, Z=1.0)), "k-",
           label="c=4, Z=1.0")
    p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
    p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
             " \Sigma_{\mathrm{HI}}$")
    p.xlim([2, 22])
    p.ylim([-4, 1])
    p.legend(loc='lower right', frameon=True)
    p.grid()
    p.tight_layout()

    save_name = "ratio_totalsigma_w_krumholzmodel_dr"
    p.savefig(paper1_figures_path("{0}_{1}pc.pdf".format(save_name,
                                                         int(dr.value))))
    p.savefig(paper1_figures_path("{0}_{1}pc.png".format(save_name,
                                                         int(dr.value))))
    p.close()

    # Gratier+16 find evidence for a dark CO component, at about ~5 Msol/pc^2.
    # Let's add this in, assuming the dark component is *only* in the CO and
    # not due to optically thick HI (some portion probably is).
    sd_dark = sd + 5 * u.solMass / u.pc**2
    sd_dark_sigma = (sd_dark * sd_sigma) / sd
    gas_ratio_dark = sd_dark.value / sd_hi.value
    gas_ratio_dark_sigma = \
        (gas_ratio_dark *
         np.sqrt((sd_dark_sigma / sd_dark)**2 +
                 (sd_sigma_hi / sd_hi)**2)).value
    log_gas_ratio_dark_sigma = gas_ratio_dark_sigma / \
        (gas_ratio_dark * np.log(10))

    total_sd_plus_dark = sd_dark.value + sd_hi.value
    total_sd_plus_dark_sigma = \
        (total_sd_plus_dark *
         np.sqrt((sd_dark_sigma / sd_dark)**2 +
                 (sd_sigma_hi / sd_hi)**2)).value

    onecolumn_figure(font_scale=1.0, fig_ratio=1.0)

    p.errorbar(total_sd_plus_dark, np.log10(gas_ratio_dark),
               yerr=log_gas_ratio_dark_sigma,
               xerr=total_sd_plus_dark_sigma, color='r', alpha=0.6, marker='o',
               label=r"H$_2$ + HI + CO-dark H$_2$")
    p.errorbar(total_sd, np.log10(gas_ratio), yerr=log_gas_ratio_sigma,
               xerr=total_sd_sigma, color='b', alpha=0.6, fmt='D',
               label=r"H$_2$ + HI")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=1, Z=0.5)), "r--",
           label="c=1, Z=0.5")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=4, Z=0.5)), "g-.",
           label="c=4, Z=0.5")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=4, Z=0.25)), "m.",
           label="c=4, Z=0.25")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=4, Z=1.0)), "k-",
           label="c=4, Z=1.0")
    p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
    p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
             " \Sigma_{\mathrm{HI}}$")
    p.xlim([2, 25])
    p.ylim([-4, 1])
    p.legend(loc='lower right', frameon=True)
    p.grid()
    p.tight_layout()

    save_name = "ratio_totalsigma_dark_w_krumholzmodel_dr"
    p.savefig(paper1_figures_path("{0}_{1}pc.pdf".format(save_name,
                                                         int(dr.value))))
    p.savefig(paper1_figures_path("{0}_{1}pc.png".format(save_name,
                                                         int(dr.value))))
    p.close()

    # But M33 has a known metallicity gradient, so we can do a bit better
    # Clumping factors should converge to 1 on 100 pc, based on the Krumholz
    # model. This isn't happening here, so let's what c needs to be for the
    # curve to intersect with each point we have.

    # Metallicity of 0.5
    clump_constz = optimize_clump_factors(total_sd, gas_ratio, Z=0.5)

    # Metallicity Gradient from Roso & Simon (2005)
    def ros_sim_metallicity(radius):
        return 10 ** (8.36 - 0.027 * radius - 8.8)

    clump_rossim = optimize_clump_factors(total_sd, gas_ratio,
                                          Z=ros_sim_metallicity(rs.value))

    # And from Bresolin 2011
    def bresolin_metallicity(radius):
        return 10 ** (8.82 - 0.03 * radius - 8.8)

    clump_bresolin = optimize_clump_factors(total_sd, gas_ratio,
                                            Z=bresolin_metallicity(rs.value))

    onecolumn_figure(font_scale=1.0)

    p.plot(rs.value[:-1], clump_constz[:-1], 'bD-', label="Z=0.5")
    p.plot(rs.value[:-1], clump_rossim[:-1], 'ro--',
           label="Rosolowsky & Simon (2005)")
    p.plot(rs.value[:-1], clump_bresolin[:-1], 'gs-.', label="Bresolin (2011)")
    p.legend(loc='best', frameon=True)
    p.ylim([-1, 10])
    p.grid()
    p.ylabel("Clumping Factor")
    p.xlabel("Radius (kpc)")
    p.tight_layout()

    save_name = "clumpfactor_krumholzmodel_dr"
    p.savefig(paper1_figures_path("{0}_{1}pc.pdf".format(save_name,
                                                         int(dr.value))))
    p.savefig(paper1_figures_path("{0}_{1}pc.png".format(save_name,
                                                         int(dr.value))))
    p.close()

    # What are the clumping factors on a per pixel basis?

    # Now load in the zeroth moments.

    hi_mom0 = fits.open(fourteenB_HI_data_path(moment0_name))[0]
    beam = Beam.from_fits_header(hi_mom0.header)
    # Convert to K km/s
    hi_mom0_data = hi_mom0.data * beam.jtok(hi_freq).value / 1000.
    hi_mom0_data = hi_mom0_data * u.K * u.km / u.s

    mom0_reproj = fits.open(iram_co21_data_path("m33.co21_iram.masked.moment0.hireprojection.fits"))
    mom0_reproj = mom0_reproj.data * u.K * u.km / u.s

    good_pts = np.where(np.isfinite(mom0_reproj))

    # Make a radius array
    radii = gal.radius(header=hi_mom0.header).to(u.kpc)
    radii_pts = radii[good_pts]

    # And the position angles
    pang = gal.position_angles(header=hi_mom0.header).to(u.deg)
    pang_pts = pang[good_pts]

    skycoord_grid = gal.skycoord_grid(header=hi_mom0.header)
    skycoord_pts = skycoord_grid[good_pts]

    # Correct for the disk inclincation
    inc = np.cos(gal.inclination)

    # 30 m beam efficiency
    beam_eff = 0.75

    # Convert the integrated intensities to surface densities.
    # hi_coldens = hi_mom0_reproj[good_pts] * hi_mass_conversion * inc
    hi_coldens = hi_mom0_data[good_pts] * hi_mass_conversion * inc

    # co_coldens = mom0[good_pts] * co21_mass_conversion * inc / beam_eff
    co_coldens = mom0_reproj[good_pts] * co21_mass_conversion * inc / beam_eff

    # Remove any NaNs in either
    nans = np.logical_or(np.isnan(hi_coldens), np.isnan(co_coldens))

    hi_coldens = hi_coldens[~nans]
    co_coldens = co_coldens[~nans]
    radii_pts = radii_pts[~nans]
    pang_pts = pang_pts[~nans]
    skycoord_pts = skycoord_pts[~nans]
    ypts = good_pts[0][~nans]
    xpts = good_pts[1][~nans]

    gas_ratio_pix = co_coldens / hi_coldens
    total_sd_pix = co_coldens + hi_coldens

    from corner import hist2d

    sds = np.arange(0.1, 60, 0.2)

    default_figure()

    hist2d(total_sd_pix.value, np.log10(gas_ratio_pix.value),
           data_kwargs={"alpha": 0.6})
    p.xlim([0, 60])
    p.ylim([-2.1, 0.7])
    p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
    p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
             " \Sigma_{\mathrm{HI}}$")
    p.errorbar(total_sd, np.log10(gas_ratio), yerr=log_gas_ratio_sigma,
               xerr=total_sd_sigma, color='b', alpha=0.6, fmt='D',
               label='100 pc radial bins')

    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=1, Z=0.5)), "r--",
           label="c=1, Z=0.5")
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=3, Z=0.5)), "g-.",
           label="c=3, Z=0.5", linewidth=4)
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=3, Z=0.5)), "m.",
           label="c=3, Z=0.5", linewidth=4)
    p.plot(sds, np.log10(krumholz_ratio_model(sds, c=2, Z=1.0)), "b-",
           label="c=2, Z=1.0", linewidth=4)
    p.legend(loc='lower right', frameon=True)

    save_name = "ratio_totalsigma_w_krumholzmodel_perpix"
    p.savefig(paper1_figures_path("{}.pdf".format(save_name)))
    p.savefig(paper1_figures_path("{}.png".format(save_name)))
    p.close()

    # Comparison with radius
    # colors = p.cm.viridis(radii_pts.value)
    # p.scatter(total_sd_pix, np.log10(gas_ratio_pix), c=colors)

    bin_width = 2
    inners = np.arange(0, 8, bin_width)
    outers = np.arange(bin_width, 8.5, bin_width)
    fig, axes = p.subplots(2, 2, sharex=True, sharey=True)
    p.subplots_adjust(hspace=0.04, wspace=0.04)

    for i, (ax, lower, upper) in enumerate(zip(axes.flatten(), inners, outers)):
        bin_pts = np.logical_and(radii_pts.value >= lower,
                                 radii_pts.value < upper)
        hist2d(total_sd_pix.value,
               np.log10(gas_ratio_pix.value),
               data_kwargs={"alpha": 0.6},
               ax=ax)
        ax.plot(total_sd_pix[bin_pts], np.log10(gas_ratio_pix[bin_pts]), "D",
                label="{0}-{1}".format(lower, upper),
                ms=3.0, mec=None, alpha=0.8)
        ax.grid()
        ax.text(40, -0.8, "{0} - {1} kpc".format(int(lower), int(upper)),
                bbox={"boxstyle": "square", "facecolor": "w"})
        ax.set_xlim([0, 60])
        ax.set_ylim([-2.1, 0.7])

        if i == 0 or i == 2:
            ax.set_ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
                          " \Sigma_{\mathrm{HI}}$")
        if i == 2 or i == 3:
            ax.set_xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")

    fig.savefig(paper1_figures_path("ratio_totalsigma_radialbins_perpix.pdf"))
    fig.savefig(paper1_figures_path("ratio_totalsigma_radialbins_perpix.png"))
    p.close()

    # p.legend()

    # Compare different properties vs. radius and PA
    fig, ax = p.subplots(1, 2, sharex=True)
    hist2d(radii_pts.value, np.log10(gas_ratio_pix.value),
           data_kwargs={"alpha": 0.6}, ax=ax[0], bins=10)
    # Make radial bins and find the median in each
    rad_bins = np.arange(0, 8.5, 0.5)
    med_ratio, bin_edges, cts = binned_statistic(radii_pts.value,
                                                 np.log10(gas_ratio_pix.value),
                                                 bins=rad_bins,
                                                 statistic=np.median)
    lower_ratio = binned_statistic(radii_pts.value,
                                   np.log10(gas_ratio_pix.value),
                                   bins=rad_bins,
                                   statistic=lambda x: np.percentile(x, 15))[0]
    upper_ratio = binned_statistic(radii_pts.value,
                                   np.log10(gas_ratio_pix.value),
                                   bins=rad_bins,
                                   statistic=lambda x: np.percentile(x, 85))[0]
    bin_cents = (bin_edges[1:] + bin_edges[:-1]) / 2.
    ax[0].errorbar(bin_cents, med_ratio, color='b', fmt='o--',
                   yerr=[med_ratio - lower_ratio, upper_ratio - med_ratio])

    ax[0].set_ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
                     " \Sigma_{\mathrm{HI}}$")
    ax[0].set_xlabel("Radius (kpc)")
    ax[0].set_xlim(0, 8.2)
    ax[1].plot(radii_pts.value, total_sd_pix, "ko",
               ms=3.0, mec=None, alpha=0.8)
    ax[1].set_ylabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
    ax[1].set_xlabel("Radius (kpc)")

    fig.savefig(paper1_figures_path("ratio_totalsigma_vs_radius_perpix.pdf"))
    fig.savefig(paper1_figures_path("ratio_totalsigma_vs_radius_perpix.png"))
    p.close()

    # fig, ax = p.subplots(1, 2)
    # # hist2d(pang_pts.value, np.log10(gas_ratio_pix.value),
    # #        data_kwargs={"alpha": 0.6}, ax=ax[0], bins=10)
    # ax[0].plot(pang_pts.value, np.log10(gas_ratio_pix.value), "ko",
    #            ms=3.0, mec=None, alpha=0.8)
    # ax[0].set_ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
    #                  " \Sigma_{\mathrm{HI}}$")
    # ax[0].set_xlabel("Position Angle (deg)")
    # ax[1].plot(pang_pts.value, total_sd_pix, "ko",
    #            ms=3.0, mec=None, alpha=0.8)
    # ax[1].set_ylabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
    # ax[1].set_xlabel("Position Angle (deg)")

    default_figure()

    # Save the lists of points in a table
    tab = Table([skycoord_pts.ra, skycoord_pts.dec, radii_pts, pang_pts,
                 hi_coldens, co_coldens, total_sd_pix, gas_ratio_pix,
                 ypts, xpts],
                names=["RA", "Dec", "Radius", "PA", "Sigma_HI", "Sigma_H2",
                       "Sigma_Total", "Ratio", "xpix", "ypix"])
    tab.write(fourteenB_HI_data_path("tables/column_densities_perpix.fits", no_check=True))
