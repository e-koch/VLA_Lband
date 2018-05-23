
'''
Use the fitted version of the Brandt curve to create a map of shear.
'''

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from cube_analysis.rotation_curves.curve_fitting import (oortA_brandt,
                                                         epifreq_brandt,
                                                         vcirc_brandt)

from paths import (fourteenB_wGBT_HI_file_dict, allfigs_path,
                   fourteenB_HI_data_wGBT_path, iram_co21_14B088_data_path)
from galaxy_params import gal_feath as gal
from constants import hi_mass_conversion, hi_freq

from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_twopanel_figure)


def shear_sd_criteria(sigma, A):
    '''
    The critical shear surface density for instabilities to occur.
    From Hunter, Elmegreen & Baker (1998) (Eq'n 2--5).
    '''
    alpha = 2.5

    crit_sd = alpha * A.to(u.km / u.s / u.kpc) * sigma.to(u.km / u.s) / \
        (np.pi * c.G)
    return crit_sd.to(u.solMass / u.pc**2)


def toomre_sd_criteria(sigma, kappa):
    '''
    Critical surface density for the gas component to be Toomre-unstable.
    Assumed a critical value at Q=1 with a correction of alpha=0.7 for finite
    disk thickness (Romeo 1992).
    '''
    alpha = 0.7

    crit_sd = alpha * kappa.to(u.km / u.s / u.kpc) * sigma.to(u.km / u.s) / \
        (np.pi * c.G)
    return crit_sd.to(u.solMass / u.pc**2)


mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0])

lwidth = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['LWidth'])[0]).to(u.km / u.s)

SurfDens_HI = mom0.quantity * mom0.beam.jtok(hi_freq) * (u.beam / u.Jy) * \
    hi_mass_conversion * np.cos(gal.inclination) / (1000. * u.m / u.km)


# Index, v_max (km/s), R_max (kpc)
pars = [0.56, 110.0, 12.0]

radii = gal.radius(header=mom0.header).to(u.kpc)

shear_map = oortA_brandt(radii.value, *pars) * u.km / u.s / u.kpc
shear_map = shear_map.to(u.km / u.s / u.pc)

# And make one of epicyclic frequency

epifreq_map = epifreq_brandt(radii.value, *pars) * u.km / u.s / u.kpc

# Make radial slice versions

radii_1D = np.arange(0., 7., 0.1) * u.kpc

shear_1D = oortA_brandt(radii_1D.value, *pars) * u.km / u.s / u.kpc
epifreq_1D = epifreq_brandt(radii_1D.value, *pars) * u.km / u.s / u.kpc

# Taking centroid stacked value.
# dispersion
sigma_gas = 6 * u.km / u.s
# Constant Sigma_HI with an exponential H2 disk
Sigma_gas = 8 * u.solMass / u.pc**2 + \
    10 * u.solMass / u.pc**2 * np.exp(-radii_1D / (2.2 * u.kpc))

Sigma_gas_map = SurfDens_HI + \
    10 * u.solMass / u.pc**2 * np.exp(-radii / (2.2 * u.kpc))

Qgas = ((sigma_gas * epifreq_1D) / (np.pi * c.G * Sigma_gas)).to(u.dimensionless_unscaled)

# Make a gas Q map assuming line widths from the second moment
# Assume that the line width from the second moment is ~50\% overestimated
Qgas_map = (lwidth * epifreq_map / (np.pi * c.G * SurfDens_HI)).to(u.dimensionless_unscaled)
Qgas_map[Qgas_map <= 0.] = np.NaN

# Exponential stellar disk. From Corbelli+14, this is about
# 242 Msol/pc^2 * exp(-R/1.7 kpc)
l_star = 1.7 * u.kpc
Sigma_star_0 = 242 * u.solMass / u.pc**2
Sigma_stellar = Sigma_star_0 * np.exp(-radii_1D / l_star)

Sigma_stellar_map = Sigma_star_0 * np.exp(- radii / l_star)

# Thin disk dispersion (corrected) from McConnachie+2006
# sigma_stellar = 12.5 * u.km / u.s

# Used in Corbelli 2003
sigma_stellar = 20 * u.km / u.s

# Equilibrium approximation used in Leroy+08
sigma_stellar_equilib = np.sqrt(2 * np.pi * c.G * l_star * Sigma_stellar / 7.3).to(u.km / u.s)

Qstar = ((sigma_stellar * epifreq_1D) / (np.pi * c.G * Sigma_stellar)).to(u.dimensionless_unscaled)
Qstar_equilib = ((sigma_stellar_equilib * epifreq_1D) / (np.pi * c.G * Sigma_stellar)).to(u.dimensionless_unscaled)

Qstar_map = ((sigma_stellar * epifreq_map) / (np.pi * c.G * Sigma_stellar_map)).to(u.dimensionless_unscaled)
Qstar_map[radii >= 10 * u.kpc] = np.NaN

prefix = (2 * sigma_gas * sigma_stellar) / (sigma_gas**2 + sigma_stellar**2)
Qeff_old = 1 / (prefix * Qstar**-1 + Qgas**-1)


def Qeff_RW11_2comp(Qs, sigma_s, Qg, sigma_g):
    '''
    Qeff from Romeo + Wiegert (2011) for two components
    '''

    prefix = (2 * sigma_g * sigma_s) / (sigma_g**2 + sigma_s**2)

    Qeff = np.zeros_like(Qs)

    Tg = 1.5
    Ts = 1.22

    Qeff[Ts * Qs >= Qg * Tg] = \
        1 / (prefix * Qs**-1 + Qg**-1)[Ts * Qs >= Qg * Tg]
    Qeff[Ts * Qs < Qg * Tg] = 1 / (prefix * Qg**-1 + Qs**-1)[Ts * Qs < Qg * Tg]

    return Qeff


Qeff = Qeff_RW11_2comp(Qstar, sigma_stellar, Qgas, sigma_gas)
Qeff_equilib = Qeff_RW11_2comp(Qstar_equilib, sigma_stellar_equilib, Qgas,
                               sigma_gas)

# Again, assume that the second moment line width is ~50\% overestimated,
# so divide by 2
Qeff_2_map = Qeff_RW11_2comp(Qstar_map, sigma_stellar, Qgas_map, lwidth)
Qeff_2_map[np.isnan(mom0.value)] = np.NaN

# Save the map
Qeff_2_hdu = fits.PrimaryHDU(Qeff_2_map.value, header=mom0.header)
Qeff_2_hdu.header['BUNIT'] = ("", "Effective Toomre Q")
Qeff_2_hdu.header['COMMENT'] = "Qeff from Romeo + Wiegert (2011). Using line"\
    " width from second moment for gas dispersion."
Qeff_2_hdu.writeto(fourteenB_HI_data_wGBT_path("M33_14B-088_HI_Qeff_map.fits",
                                               no_check=True))

onecolumn_figure()
plt.axhline(1, color='k', linestyle='--', alpha=0.5)
plt.plot(radii_1D, Qgas, label=r'$Q_{\rm gas}$')
plt.plot(radii_1D, Qstar, label=r'$Q_{\rm star}$')
plt.plot(radii_1D, Qstar_equilib, label=r'$Q_{\rm star}$ Equilib.')
plt.plot(radii_1D, Qeff, label=r'$Q_{\rm eff}$')
plt.plot(radii_1D, Qeff_old, label=r'$Q_{\rm eff}$ OLD')
plt.plot(radii_1D, Qeff_equilib, label=r'$Q_{\rm eff}$ Equilib.')
plt.grid()
plt.legend(frameon=True, loc='upper right')
plt.ylabel(r"Toomre $Q$")
plt.xlabel("R (kpc)")
plt.ylim([0, 4])
plt.tight_layout()

plt.savefig(allfigs_path("toomre_q_radius.png"))
plt.savefig(allfigs_path("toomre_q_radius.pdf"))

plt.close()


# Now following Romeo + Mogotsi 2017, treat the HI and CO separately.
# I'll load in the centroid-centered profile values here.

def Qeff_RW11_3comp(Qs, sigma_s, QHI, sigma_HI, QH2, sigma_H2):
    '''
    Qeff from Romeo + Wiegert (2011) for two components
    '''

    Qeff = np.zeros_like(Qs).value

    Tg = 1.5
    Ts = 1.22

    for i in range(len(Qs)):

        try:
            sigma_s_i = sigma_s[i]
        except TypeError:
            sigma_s_i = sigma_s

        try:
            sigma_HI_i = sigma_HI[i]
        except TypeError:
            sigma_HI_i = sigma_HI

        try:
            sigma_H2_i = sigma_H2[i]
        except TypeError:
            sigma_H2_i = sigma_H2

        # Find the min of T_i Q_i
        sigmas = np.array([sigma_s_i.value, sigma_HI_i.value,
                           sigma_H2_i.value])

        weights = np.array([Ts * Qs[i].value, Tg * QHI[i].value,
                            Tg * QH2[i].value])

        sigma_m = sigmas[np.argmin(weights)] * u.km / u.s

        W_s = (2 * sigma_m * sigma_s_i) / (sigma_m**2 + sigma_s_i**2)
        W_hi = (2 * sigma_m * sigma_HI_i) / (sigma_m**2 + sigma_HI_i**2)
        W_h2 = (2 * sigma_m * sigma_H2_i) / (sigma_m**2 + sigma_H2_i**2)

        term1 = W_s / (Ts * Qs[i])
        term2 = W_hi / (Tg * QHI[i])
        term3 = W_h2 / (Tg * QH2[i])

        Qeff[i] = (1 / (term1 + term2 + term3)).value

    return Qeff


hi_tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_radial_500pc.csv"))
co_tab = Table.read(iram_co21_14B088_data_path("tables/co_hwhm_totalprof_fits_radial_500pc.csv"))
bin_centers = (np.arange(0., 6.6, 0.5) + 0.25) * u.kpc

Sigma_HI = 8 * u.solMass / u.pc**2
Sigma_H2 = 10 * u.solMass / u.pc**2 * np.exp(-bin_centers / (2.2 * u.kpc))

sigma_hi = hi_tab['centsub_sigma'] * u.km / u.s
sigma_h2 = co_tab['centsub_sigma'] * u.km / u.s

epifreq_bin = epifreq_brandt(bin_centers.value, *pars) * u.km / u.s / u.kpc
shear_bin = oortA_brandt(bin_centers.value, *pars) * u.km / u.s / u.kpc

QHI = ((sigma_hi * epifreq_bin) / (np.pi * c.G * Sigma_HI)).to(u.dimensionless_unscaled)
QH2 = ((sigma_h2 * epifreq_bin) / (np.pi * c.G * Sigma_H2)).to(u.dimensionless_unscaled)

Sigma_stellar_bin = Sigma_star_0 * np.exp(-bin_centers / l_star)
sigma_stellar_equilib_bin = np.sqrt(2 * np.pi * c.G * l_star * Sigma_stellar_bin / 7.3).to(u.km / u.s)

Qstar_bin = ((sigma_stellar * epifreq_bin) / (np.pi * c.G * Sigma_stellar_bin)).to(u.dimensionless_unscaled)
Qstar_bin_equilib = ((sigma_stellar_equilib_bin * epifreq_bin) / (np.pi * c.G * Sigma_stellar_bin)).to(u.dimensionless_unscaled)

Qeff_3 = Qeff_RW11_3comp(Qstar_bin, sigma_stellar, QHI, sigma_hi, QH2, sigma_h2)
Qeff_3_equilib = Qeff_RW11_3comp(Qstar_bin_equilib, sigma_stellar_equilib_bin,
                                 QHI, sigma_hi, QH2, sigma_h2)

onecolumn_figure()
plt.axhline(1, color='k', linestyle='--', alpha=0.5)
plt.plot(bin_centers, QHI, label=r'$Q_{\rm HI}$', drawstyle='steps-mid')
plt.plot(bin_centers, QH2, label=r'$Q_{\rm H2}$', drawstyle='steps-mid')
plt.plot(bin_centers, Qstar_bin, label=r'$Q_{\rm star}$',
         drawstyle='steps-mid')
plt.plot(bin_centers, Qstar_bin_equilib, label=r'$Q_{\rm star}$ Equilib.',
         drawstyle='steps-mid')
plt.plot(bin_centers, Qeff_3, label=r'$Q_{\rm eff}$',
         drawstyle='steps-mid')
plt.plot(bin_centers, Qeff_3_equilib, label=r'$Q_{\rm eff}$ Equilib.',
         drawstyle='steps-mid')
plt.grid()
plt.legend(frameon=True)
plt.ylim([0, 8])
plt.ylabel(r"$Q$")
plt.xlabel("R (kpc)")
plt.tight_layout()

plt.savefig(allfigs_path("toomre_q3_bins.png"))
plt.savefig(allfigs_path("toomre_q3_bins.pdf"))

plt.ylim([-0.2, 2])
plt.legend(frameon=True, ncol=2)

plt.savefig(allfigs_path("toomre_q3_bins_zoom.png"))
plt.savefig(allfigs_path("toomre_q3_bins_zoom.pdf"))

# Make a full map of Qeff
Qeff_map = np.zeros_like(radii.value)

bin_low = np.arange(0., 6.6, 0.5) * u.kpc
bin_high = np.arange(0.5, 7.1, 0.5) * u.kpc

for idx, (low, high) in enumerate(zip(bin_low, bin_high)):

    rad_mask = np.logical_and(radii >= low, radii <= high)

    Qeff_map[rad_mask] = Qeff_3[idx]

Qeff_map[Qeff_map == 0.] = np.NaN

# Compare surface density thresholds

onecolumn_figure()
plt.plot(radii_1D, Sigma_gas, label='Gas')
plt.plot(radii_1D, shear_sd_criteria(sigma_gas, shear_1D), label='Shear')
plt.plot(radii_1D, toomre_sd_criteria(sigma_gas, epifreq_1D), label='Toomre')
plt.grid()
plt.legend(frameon=True, loc='upper right')
plt.ylabel(r"$\Sigma_{\rm gas}$ (M$_\odot$ pc$^{-2}$)")
plt.xlabel("R (kpc)")
plt.tight_layout()

plt.savefig(allfigs_path("gas_crit_sd_radius.png"))
plt.savefig(allfigs_path("gas_crit_sd_radius.pdf"))

plt.close()

# And for the binned values
twocolumn_twopanel_figure()
plt.subplot(131)
plt.plot(bin_centers, [Sigma_HI.value] * 14, label='HI',
         drawstyle='steps-mid')
plt.plot(bin_centers, shear_sd_criteria(sigma_hi, shear_bin), label='Shear',
         drawstyle='steps-mid')
plt.plot(bin_centers,
         toomre_sd_criteria(sigma_h2, epifreq_bin), label='Toomre',
         drawstyle='steps-mid')
plt.grid()
plt.legend(frameon=True, loc='upper right')
plt.ylabel(r"$\Sigma_{\rm gas}$ (M$_\odot$ pc$^{-2}$)")
plt.xlabel("R (kpc)")
plt.subplot(132)
plt.plot(bin_centers, Sigma_H2, label='H2',
         drawstyle='steps-mid')
plt.plot(bin_centers, shear_sd_criteria(sigma_h2, shear_bin), label='Shear',
         drawstyle='steps-mid')
plt.plot(bin_centers,
         toomre_sd_criteria(sigma_h2, epifreq_bin), label='Toomre',
         drawstyle='steps-mid')
plt.grid()
plt.legend(frameon=True, loc='upper right')
plt.ylabel(r"$\Sigma_{\rm gas}$ (M$_\odot$ pc$^{-2}$)")
plt.xlabel("R (kpc)")
plt.subplot(133)
plt.plot(bin_centers, Sigma_H2 + Sigma_HI, label='H2',
         drawstyle='steps-mid')
plt.plot(bin_centers, shear_sd_criteria(sigma_hi, shear_bin), label='Shear',
         drawstyle='steps-mid')
plt.plot(bin_centers,
         toomre_sd_criteria(sigma_hi, epifreq_bin), label='Toomre',
         drawstyle='steps-mid')
plt.grid()
plt.legend(frameon=True, loc='upper right')
plt.ylabel(r"$\Sigma_{\rm gas}$ (M$_\odot$ pc$^{-2}$)")
plt.xlabel("R (kpc)")
plt.tight_layout()

# Location of LRs
Omega = vcirc_brandt(radii_1D.value, *pars) * u.km / u.s / radii_1D

onecolumn_figure()
plt.plot(radii_1D, Omega, label=r'$\Omega$')
plt.plot(radii_1D, Omega - epifreq_1D / 2., label=r'$\Omega - \kappa / 2$')
plt.plot(radii_1D, Omega + epifreq_1D / 2., label=r'$\Omega + \kappa / 2$')
plt.grid()
plt.legend(frameon=True, loc='upper right')
plt.ylabel(r"$\Sigma_{\rm gas}$ (M$_\odot$ pc$^{-2}$)")
plt.xlabel("R (kpc)")
plt.tight_layout()
