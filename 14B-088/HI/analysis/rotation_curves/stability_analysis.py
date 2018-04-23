
'''
Use the fitted version of the Brandt curve to create a map of shear.
'''

from astropy.io import fits
import astropy.units as u
import astropy.constants as c
import numpy as np
import matplotlib.pyplot as plt

from cube_analysis.rotation_curves.curve_fitting import (oortA_brandt,
                                                         epifreq_brandt,
                                                         vcirc_brandt)

from paths import fourteenB_wGBT_HI_file_dict, allfigs_path
from galaxy_params import gal_feath as gal

from plotting_styles import default_figure, onecolumn_figure


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
    Assumed a critical value at Q=1 (this is alpha).
    '''
    alpha = 1.0

    crit_sd = alpha * kappa.to(u.km / u.s / u.kpc) * sigma.to(u.km / u.s) / \
        (np.pi * c.G)
    return crit_sd.to(u.solMass / u.pc**2)


mom0 = fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0]

# Index, v_max (km/s), R_max (kpc)
pars = [0.56, 110.0, 12.0]

radii = gal.radius(header=mom0.header).to(u.kpc)

shear_map = oortA_brandt(radii.value, *pars) * u.km / u.s / u.kpc
shear_map = shear_map.to(u.km / u.s / u.pc)

# And make one of epicyclic frequency

epifreq_map = epifreq_brandt(radii.value, *pars) * u.km / u.s

# Make radial slice versions

radii_1D = np.arange(0., 7., 0.1) * u.kpc

shear_1D = oortA_brandt(radii_1D.value, *pars) * u.km / u.s / u.kpc
epifreq_1D = epifreq_brandt(radii_1D.value, *pars) * u.km / u.s / u.kpc


# Taking centroid stacked value.
# dispersion
sigma_gas = 10 * u.km / u.s
# Constant Sigma_HI with an exponential H2 disk
Sigma_gas = 10 * u.solMass / u.pc**2 + \
    10 * u.solMass / u.pc**2 * np.exp(-radii_1D / (2.2 * u.kpc))

Qgas = ((sigma_gas * epifreq_1D) / (np.pi * c.G * Sigma_gas)).to(u.dimensionless_unscaled)

# Exponential stellar disk. From Corbelli+14, this is about
# 200 Msol/pc^2 * exp(-R/2 kpc)
Sigma_stellar = 200 * u.solMass / u.pc**2 * np.exp(-radii_1D / (2 * u.kpc))

# Thin disk dispersion from McConnachie+2006
sigma_stellar = 16 * u.km / u.s

Qstar = ((sigma_stellar * epifreq_1D) / (np.pi * c.G * Sigma_stellar)).to(u.dimensionless_unscaled)

prefix = (2 * sigma_gas * sigma_stellar) / (sigma_gas**2 + sigma_stellar**2)
Qeff = 1 / (prefix * Qstar**-1 + Qgas**-1)

onecolumn_figure()
plt.axhline(1, color='k', linestyle='--', alpha=0.5)
plt.plot(radii_1D, Qgas, label=r'$Q_{\rm gas}$')
plt.plot(radii_1D, Qstar, label=r'$Q_{\rm star}$')
plt.plot(radii_1D, Qeff, label=r'$Q_{\rm eff}$')
plt.grid()
plt.legend(frameon=True, loc='upper right')
plt.ylabel(r"Toomre $Q$")
plt.xlabel("R (kpc)")
plt.ylim([0, 4])
plt.tight_layout()

plt.savefig(allfigs_path("toomre_q_radius.png"))
plt.savefig(allfigs_path("toomre_q_radius.pdf"))

plt.close()

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
