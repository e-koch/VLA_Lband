
from astropy.io import fits
from galaxies import Galaxy
import matplotlib.pyplot as p
import os
from astropy.coordinates import Angle
from astropy import units as u
import numpy as np
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from astropy.wcs import WCS
from astropy.table import Table

# Import from above.
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

from HI_radial_profile import surfdens_radial_profile
from paths import (iram_co21_data_path, fourteenB_HI_data_path,
                   paper1_figures_path, c_hi_analysispath)
from constants import moment0_name, cube_name

from krumholz_2009_model import krumholz_ratio_model

'''
Create the surface density profile in CO(2-1), assuming a factor to convert
to the H2 mass.
'''

co21_mass_conversion = 6.7 * (u.M_sun / u.pc ** 2) / (u.K * u.km / u.s)

# IRAM beam efficiency
beam_eff = 0.75

# Set the radial disk widths to bin over
# dr = 500 * u.pc
dr = 100 * u.pc
# dr = 300 * u.pc

# Load the moment 0
cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
cube = cube.with_mask(cube > 0.1 * u.K)

direc_hi = "/home/eric/MyRAID/M33/14B-088/HI/full_imaging/"
mom0_hi = fits.open(fourteenB_HI_data_path(moment0_name))[0]
hi_cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

g = Galaxy("M33")
radii = g.radius(header=cube.header)
# Edge effects are really awful in this map. Ignore the edges by masking
# beyond 6 kpc. This is really close to the edge of the data anyways, and
# honestly results beyond this point shouldn't be trusted...
cube = cube.with_mask(radii < 6. * u.kpc)

# mom0 = fits.open(os.path.join(direc, "m33.ico.fits"))[0]

# mom0_data = mom0.data.squeeze() * (mom0.data.squeeze() > 1.0) * u.K
mom0 = cube.moment0()

rs, sd, sd_sigma = surfdens_radial_profile(g, cube=cube, mom0=mom0,
                                           max_rad=6 * u.kpc, dr=dr)
# Correct for beam efficiency
sd /= beam_eff
sd_sigma /= beam_eff

rs_n, sd_n, sd_sigma_n = \
    surfdens_radial_profile(g, cube=cube, mom0=mom0,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                             -0.5 * np.pi * u.rad]),
                            max_rad=6 * u.kpc, dr=dr)
sd_n /= beam_eff
sd_sigma_n /= beam_eff

rs_s, sd_s, sd_sigma_s = \
    surfdens_radial_profile(g, cube=cube, mom0=mom0,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            max_rad=6 * u.kpc, dr=dr)
sd_n /= beam_eff
sd_sigma_n /= beam_eff

p.ioff()
p.errorbar(rs.value, np.log10(sd.value),
           yerr=0.434 * sd_sigma.value / sd.value, fmt="-", color="b",
           label=r"H$_2$", drawstyle='steps-mid')
p.ylabel(r" log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"R (kpc)")
# p.legend(loc='best')
p.grid("on")

p.savefig(paper1_figures_path("M33_Sigma_profile_co21_dr_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("M33_Sigma_profile_co21_dr_{}pc.png".format(dr.value)))
p.close()
# p.show()

# Show the north vs south profiles
p.plot(rs.value, np.log10(sd.value), "k-.", drawstyle='steps-mid',
       label="Total")
p.errorbar(rs_n.value, np.log10(sd_n.value),
           yerr=0.434 * sd_sigma_n.value / sd_n.value, fmt="-", color="b",
           label="North", drawstyle='steps-mid')
p.errorbar(rs_s.value, np.log10(sd_s.value),
           yerr=0.434 * sd_sigma_s.value / sd_s.value, fmt="--", color="g",
           label="South", drawstyle='steps-mid')

# p.plot(rs_n.value, sd_n.value, "bD-", label="North")
# p.plot(rs_s.value, sd_s.value, "go-", label="South")
p.ylabel(r"log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"R (kpc)")
p.legend(loc='best')
p.grid("on")

p.savefig(paper1_figures_path("M33_Sigma_profile_co21_N_S_dr_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("M33_Sigma_profile_co21_N_S_dr_{}pc.png".format(dr.value)))
p.close()

# p.show()

# Now get the HI profile on the same scales
proj = Projection(mom0_hi.data * u.Jy * u.m / u.s, meta={'beam': hi_cube.beam},
                  wcs=WCS(cube[0].header))
rs_hi, sd_hi, sd_sigma_hi = surfdens_radial_profile(g, cube=hi_cube,
                                                    mom0=proj,
                                                    max_rad=6 * u.kpc, dr=dr)
# Apply scaling factor
# sd_hi /= 1.45
# sd_sigma_hi /= 1.45

# Overplot these two.
p.errorbar(rs.value, np.log10(sd.value),
           yerr=0.434 * sd_sigma.value / sd.value, fmt="-", color="b",
           label=r"H$_2$", drawstyle='steps-mid')
p.errorbar(rs_hi.value, np.log10(sd_hi.value),
           yerr=0.434 * sd_sigma_hi.value / sd_hi.value, fmt="--", color="g",
           label=r"HI", drawstyle='steps-mid')
p.ylabel(r" log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"R (kpc)")
p.legend(loc='best')
p.grid("on")

p.savefig(paper1_figures_path("M33_Sigma_profile_hi_co21_dr_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("M33_Sigma_profile_hi_co21_dr_{}pc.png".format(dr.value)))
p.close()

# p.show()

# Now plot their ratio against the total gas surface density
gas_ratio = sd.value / sd_hi.value
total_sd = sd.value + sd_hi.value
total_sd_sigma = total_sd * np.sqrt((sd_sigma / sd) ** 2 + (sd_sigma_hi / sd_hi)).value

# Overplot the Krumholz model with a few different clumping factors.
# Theoretically, c -> 1 at a resolution of 100 pc. but I'm finding a better
# match when c=4-6. The model is supposed to take metallicity into account,
# but maybe the gradient is causing some issues? Schruba+11 finds c~2 for their
# entire sample, with a lot of scatter
sds = np.arange(1, 40, 0.2)

p.semilogy(total_sd, gas_ratio, 'bD')
p.plot(sds, krumholz_ratio_model(sds, c=2, Z=0.5), "r--", label="c=2, Z=0.5")
p.plot(sds, krumholz_ratio_model(sds, c=4, Z=0.5), "g-.", label="c=4, Z=0.5")
p.plot(sds, krumholz_ratio_model(sds, c=4, Z=0.25), "m.", label="c=4, Z=0.25")
p.plot(sds, krumholz_ratio_model(sds, c=4, Z=1.0), "k-", label="c=4, Z=1.0")
p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} / \Sigma_{\mathrm{HI}}$")
p.xlim([2, 22])
p.ylim([1e-4, 10])
p.legend(loc='upper left')
p.grid()

p.savefig(paper1_figures_path("ratio_totalsigma_w_krumholzmodel_dr_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("ratio_totalsigma_w_krumholzmodel_dr_{}pc.png".format(dr.value)))
p.close()

# p.show()

# Also plot the total gas surface density against the stellar surface density
# from Corbelli
corbelli = Table.read(c_hi_analysispath("rotation_curves/corbelli_rotation_curve.csv"))

p.semilogy(rs.value, total_sd,
           linestyle="-", color="b",
           label="Gas", drawstyle='steps-mid')
p.semilogy(corbelli["R"][corbelli["R"] <= 6.5],
           corbelli["SigmaStellar"][corbelli["R"] <= 6.5], "g--",
           drawstyle='steps-mid',
           label="Stars")
p.ylabel(r"log $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"R (kpc)")
p.legend(loc='best')
p.grid()

p.savefig(paper1_figures_path("M33_Sigma_profile_gas_stars_corbelli_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("M33_Sigma_profile_gas_stars_corbelli_{}pc.png".format(dr.value)))
p.close()


# p.show()

# Finally, let's calculate some clumping factors a la Leroy+13
rs_m, sd_m, sd_sigma_m = surfdens_radial_profile(g, cube=cube, mom0=mom0,
                                                 max_rad=6 * u.kpc, dr=dr,
                                                 weight_type='mass')
# Correct for beam efficiency
sd_m /= beam_eff
sd_sigma_m /= beam_eff

rs_hi_m, sd_hi_m, sd_sigma_hi_m = \
    surfdens_radial_profile(g, cube=hi_cube,
                            mom0=proj,
                            max_rad=6 * u.kpc, dr=dr,
                            weight_type='mass')
# Apply scaling factor
# sd_hi_m /= 1.45
# sd_sigma_hi_m /= 1.45

#
p.errorbar(np.log10(sd.value), np.log10(sd_m.value),
           xerr=0.434 * sd_sigma.value / sd.value,
           yerr=0.434 * sd_sigma_m.value / sd_m.value,
           fmt="o", color="g", label="H$_2$")
p.errorbar(np.log10(sd_hi.value), np.log10(sd_hi_m.value),
           xerr=0.434 * sd_sigma_hi.value / sd_hi.value,
           yerr=0.434 * sd_sigma_hi_m.value / sd_hi_m.value,
           fmt="D", color="b", label="HI")
equality = np.arange(-2.5, 2, 0.1)
p.plot(equality, equality, 'k--')
p.ylabel(r"log Mass-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"log Area-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
p.ylim([0.35, 1.9])
p.xlim([-2.1, 1.4])
p.legend(loc='upper left')
p.show()
# The H2 (ie CO) is all over the place, and HI is clustered together and hard to see.
# Make an HI only
p.errorbar(np.log10(sd_hi.value), np.log10(sd_hi_m.value),
           xerr=0.434 * sd_sigma_hi.value / sd_hi.value,
           yerr=0.434 * sd_sigma_hi_m.value / sd_hi_m.value,
           fmt="D", color="b", label="HI")
equality = np.arange(-2.5, 2, 0.1)
p.plot(equality, equality, 'k--')
p.ylabel(r"log Mass-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"log Area-Weighted $\Sigma$ / (M$_{\odot}$ pc$^{-2}$)")
p.ylim([0.72, 1.02])
p.xlim([0.72, 1.02])
# p.legend(loc='upper left')

p.savefig(paper1_figures_path("area_weighted_vs_mass_weighted_dr_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("area_weighted_vs_mass_weighted_dr_{}pc.png".format(dr.value)))
p.close()

# p.show()


clump_co = sd_m / sd
clump_hi = sd_hi_m / sd_hi

p.ion()
