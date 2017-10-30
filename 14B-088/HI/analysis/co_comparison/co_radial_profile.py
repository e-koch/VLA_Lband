
from astropy.io import fits
import matplotlib.pyplot as p
from astropy.coordinates import Angle
from astropy import units as u
import numpy as np
from spectral_cube.lower_dimensional_structures import Projection
from astropy.table import Table
import os
from os.path import exists
from os.path import join as osjoin

from cube_analysis.profiles import surfdens_radial_profile

from paths import (iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path,
                   allfigs_path, c_hi_analysispath,
                   fourteenB_wGBT_HI_file_dict)
from constants import co21_mass_conversion, hi_mass_conversion
from galaxy_params import gal_feath as gal
from plotting_styles import onecolumn_figure, default_figure

'''
Create the surface density profile in CO(2-1), assuming a factor to convert
to the H2 mass.
'''

fig_path = osjoin(allfigs_path(""), "co_vs_hi")
if not exists(fig_path):
    os.mkdir(fig_path)

fig_co_path = osjoin(allfigs_path(""), "CO21_properties")
if not exists(fig_co_path):
    os.mkdir(fig_co_path)

# IRAM beam efficiency
beam_eff = 0.75

# Set the radial disk widths to bin over
dr = 100 * u.pc

# Load the moment 0
co_mom0 = Projection.from_hdu(fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits"))[0])
co_mom0 = co_mom0.to(u.K * u.km / u.s) / beam_eff

hi_mom0 = \
    Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0])

rs, sd, sd_sigma = surfdens_radial_profile(gal, mom0=co_mom0,
                                           max_rad=7 * u.kpc, dr=dr,
                                           mass_conversion=co21_mass_conversion)

rs_n, sd_n, sd_sigma_n = \
    surfdens_radial_profile(gal, mom0=co_mom0,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                             -0.5 * np.pi * u.rad]),
                            max_rad=7 * u.kpc, dr=dr,
                            mass_conversion=co21_mass_conversion)

rs_s, sd_s, sd_sigma_s = \
    surfdens_radial_profile(gal, mom0=co_mom0,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            max_rad=7 * u.kpc, dr=dr,
                            mass_conversion=co21_mass_conversion)

onecolumn_figure()

p.errorbar(rs.value, np.log10(sd.value),
           yerr=0.434 * sd_sigma.value / sd.value, fmt="-",
           label=r"H$_2$", drawstyle='steps-mid')
p.ylabel(r" log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.grid("on")

p.tight_layout()

p.savefig(osjoin(fig_co_path, "M33_Sigma_profile_co21_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(osjoin(fig_co_path, "M33_Sigma_profile_co21_dr_{}pc.png".format(int(dr.value))))
p.close()

# Show the north vs south profiles
p.plot(rs.value, np.log10(sd.value), "-.", drawstyle='steps-mid',
       label="Total")
p.errorbar(rs_n.value, np.log10(sd_n.value),
           yerr=0.434 * sd_sigma_n.value / sd_n.value, fmt="-",
           label="North", drawstyle='steps-mid')
p.errorbar(rs_s.value, np.log10(sd_s.value),
           yerr=0.434 * sd_sigma_s.value / sd_s.value, fmt="--",
           label="South", drawstyle='steps-mid')

# p.plot(rs_n.value, sd_n.value, "bD-", label="North")
# p.plot(rs_s.value, sd_s.value, "go-", label="South")
p.ylabel(r"log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid("on")

p.savefig(osjoin(fig_co_path, "M33_Sigma_profile_co21_N_S_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(osjoin(fig_co_path, "M33_Sigma_profile_co21_N_S_dr_{}pc.png".format(int(dr.value))))
p.close()

# p.show()

# Now get the HI profile on the same scales
rs_hi, sd_hi, sd_sigma_hi = \
    surfdens_radial_profile(gal,
                            mom0=hi_mom0,
                            max_rad=7 * u.kpc, dr=dr,
                            mass_conversion=hi_mass_conversion)

# Overplot these two.
p.errorbar(rs.value, np.log10(sd.value),
           yerr=0.434 * sd_sigma.value / sd.value, fmt="-",
           label=r"H$_2$", drawstyle='steps-mid')
p.errorbar(rs_hi.value, np.log10(sd_hi.value),
           yerr=0.434 * sd_sigma_hi.value / sd_hi.value, fmt="--",
           label=r"HI", drawstyle='steps-mid')
p.ylabel(r" log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid("on")
p.tight_layout()

p.savefig(osjoin(fig_path, "M33_Sigma_profile_hi_co21_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(osjoin(fig_path, "M33_Sigma_profile_hi_co21_dr_{}pc.png".format(int(dr.value))))
p.close()

# Save the radial profiles in a FITS table
table = Table([rs, sd, sd_sigma, sd_hi, sd_sigma_hi],
              names=('Radius', "CO_Sigma", "CO_Sigma_std", "HI_Sigma",
                     "HI_Sigma_std"))
table.write(fourteenB_HI_data_wGBT_path("tables/co21_hi_radialprofiles_{}pc.fits".format(int(dr.value)),
                                        no_check=True),
            overwrite=True)



# Also plot the total gas surface density against the stellar surface density
# from Corbelli
corbelli = Table.read(c_hi_analysispath("rotation_curves/corbelli_rotation_curve.csv"))

p.semilogy(rs.value, (sd_hi + sd).value,
           linestyle="-",
           label="Gas", drawstyle='steps-mid')
p.semilogy(corbelli["R"][corbelli["R"] <= 6.5],
           corbelli["SigmaStellar"][corbelli["R"] <= 6.5], "g--",
           drawstyle='steps-mid',
           label="Stars")
p.ylabel(r"log $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid()
p.tight_layout()

p.savefig(osjoin(fig_path, "M33_Sigma_profile_gas_stars_corbelli_{}pc.pdf".format(int(dr.value))))
p.savefig(osjoin(fig_path, "M33_Sigma_profile_gas_stars_corbelli_{}pc.png".format(int(dr.value))))
p.close()

# Finally, let's calculate some clumping factors as in Leroy+13
# rs_m, sd_m, sd_sigma_m = surfdens_radial_profile(gal, mom0=co_mom0,
#                                                  max_rad=7 * u.kpc, dr=dr,
#                                                  weight_type='mass',
#                                                  mass_conversion=co21_mass_conversion)

# rs_hi_m, sd_hi_m, sd_sigma_hi_m = \
#     surfdens_radial_profile(gal,
#                             mom0=hi_mom0,
#                             max_rad=7 * u.kpc, dr=dr,
#                             weight_type='mass')

# p.errorbar(np.log10(sd.value), np.log10(sd_m.value),
#            xerr=0.434 * sd_sigma.value / sd.value,
#            yerr=0.434 * sd_sigma_m.value / sd_m.value,
#            fmt="o", label="H$_2$")
# p.errorbar(np.log10(sd_hi.value), np.log10(sd_hi_m.value),
#            xerr=0.434 * sd_sigma_hi.value / sd_hi.value,
#            yerr=0.434 * sd_sigma_hi_m.value / sd_hi_m.value,
#            fmt="D", label="HI")
# equality = np.arange(-2.5, 2, 0.1)
# p.plot(equality, equality, 'k--')
# p.ylabel(r"log Mass-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
# p.xlabel(r"log Area-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
# p.ylim([0.25, 1.9])
# p.xlim([-2.1, 1.4])
# p.legend(loc='upper left', frameon=True)
# p.savefig(allfigs_path(osjoin(fig_path, "hi_co_area_weighted_vs_mass_weighted_dr_{}pc.pdf".format(int(dr.value)))))
# p.savefig(allfigs_path(osjoin(fig_path, "hi_co_area_weighted_vs_mass_weighted_dr_{}pc.png".format(int(dr.value)))))
# p.tight_layout()
# p.close()
# The H2 (ie CO) is all over the place, and HI is clustered together and hard
# to see. Make an HI only
# p.errorbar(np.log10(sd_hi.value), np.log10(sd_hi_m.value),
#            xerr=0.434 * sd_sigma_hi.value / sd_hi.value,
#            yerr=0.434 * sd_sigma_hi_m.value / sd_hi_m.value,
#            fmt="D", label="HI")
# equality = np.arange(-2.5, 2, 0.1)
# p.plot(equality, equality, 'k--')
# p.ylabel(r"log Mass-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
# p.xlabel(r"log Area-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
# p.ylim([0.65, 1.0])
# p.xlim([0.65, 0.9])
# p.tight_layout()

# p.savefig(allfigs_path(osjoin(fig_path, "area_weighted_vs_mass_weighted_dr_{}pc.pdf".format(int(dr.value)))))
# p.savefig(allfigs_path(osjoin(fig_path, "area_weighted_vs_mass_weighted_dr_{}pc.png".format(int(dr.value)))))
# p.close()

# clump_co = sd_m / sd
# clump_hi = sd_hi_m / sd_hi

default_figure()
