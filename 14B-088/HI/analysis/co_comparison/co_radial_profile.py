
from astropy.io import fits
import matplotlib.pyplot as p
from astropy.coordinates import Angle
from astropy import units as u
import numpy as np
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.cube_utils import average_beams
from astropy.wcs import WCS
from astropy.table import Table

from cube_analysis.profiles import surfdens_radial_profile

from paths import (iram_co21_data_path, fourteenB_HI_data_path,
                   paper1_figures_path, c_hi_analysispath)
from constants import (moment0_name, cube_name, mask_name,
                       co21_mass_conversion,
                       hi_mass_conversion)
from galaxy_params import gal
from plotting_styles import onecolumn_figure, default_figure

'''
Create the surface density profile in CO(2-1), assuming a factor to convert
to the H2 mass.
'''

default_figure()

# IRAM beam efficiency
beam_eff = 0.75

# Set the radial disk widths to bin over
# dr = 500 * u.pc
dr = 100 * u.pc
# dr = 300 * u.pc

# Load the moment 0
cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
# cube = cube.with_mask(cube > 0.1 * u.K)

mom0_hi = fits.open(fourteenB_HI_data_path(moment0_name))[0]
hi_cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))
mask = fits.open(fourteenB_HI_data_path(mask_name))[0]
hi_cube = hi_cube.with_mask(mask.data > 0)

radii = gal.radius(header=cube.header)
# Edge effects are really awful in this map. Ignore the edges by masking
# beyond 6 kpc. This is really close to the edge of the data anyways, and
# honestly results beyond this point shouldn't be trusted...
cube = cube.with_mask(radii < 6. * u.kpc)

# mom0 = fits.open(os.path.join(direc, "m33.ico.fits"))[0]

# mom0_data = mom0.data.squeeze() * (mom0.data.squeeze() > 1.0) * u.K
mom0 = cube.moment0().to(u.K * u.km / u.s)

rs, sd, sd_sigma = surfdens_radial_profile(gal, cube=cube, mom0=mom0,
                                           max_rad=6 * u.kpc, dr=dr,
                                           mass_conversion=co21_mass_conversion)
# Correct for beam efficiency
sd /= beam_eff
sd_sigma /= beam_eff

rs_n, sd_n, sd_sigma_n = \
    surfdens_radial_profile(gal, cube=cube, mom0=mom0,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                             -0.5 * np.pi * u.rad]),
                            max_rad=6 * u.kpc, dr=dr,
                            mass_conversion=co21_mass_conversion)

sd_n /= beam_eff
sd_sigma_n /= beam_eff

rs_s, sd_s, sd_sigma_s = \
    surfdens_radial_profile(gal, cube=cube, mom0=mom0,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            max_rad=6 * u.kpc, dr=dr,
                            mass_conversion=co21_mass_conversion)

sd_s /= beam_eff
sd_sigma_s /= beam_eff

p.errorbar(rs.value, np.log10(sd.value),
           yerr=0.434 * sd_sigma.value / sd.value, fmt="-", color="b",
           label=r"H$_2$", drawstyle='steps-mid')
p.ylabel(r" log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
# p.legend(loc='best')
p.grid("on")

p.savefig(paper1_figures_path("M33_Sigma_profile_co21_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(paper1_figures_path("M33_Sigma_profile_co21_dr_{}pc.png".format(int(dr.value))))
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
p.xlabel(r"Radius (kpc)")
p.legend(loc='best')
p.grid("on")

p.savefig(paper1_figures_path("M33_Sigma_profile_co21_N_S_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(paper1_figures_path("M33_Sigma_profile_co21_N_S_dr_{}pc.png".format(int(dr.value))))
p.close()

# p.show()

# Now get the HI profile on the same scales
proj = Projection(mom0_hi.data * u.Jy * u.m / u.s, meta={'beam': average_beams(hi_cube.beams)},
                  wcs=WCS(cube[0].header))
rs_hi, sd_hi, sd_sigma_hi = surfdens_radial_profile(gal, cube=hi_cube,
                                                    mom0=proj,
                                                    max_rad=6 * u.kpc, dr=dr,
                                                    beam=average_beams(hi_cube.beams),
                                                    mass_conversion=hi_mass_conversion)
# Apply scaling factor
# sd_hi /= 1.45
# sd_sigma_hi /= 1.45

# Overplot these two.
onecolumn_figure(font_scale=1.0)

p.errorbar(rs.value, np.log10(sd.value),
           yerr=0.434 * sd_sigma.value / sd.value, fmt="-", color="b",
           label=r"H$_2$", drawstyle='steps-mid')
p.errorbar(rs_hi.value, np.log10(sd_hi.value),
           yerr=0.434 * sd_sigma_hi.value / sd_hi.value, fmt="--", color="g",
           label=r"HI", drawstyle='steps-mid')
p.ylabel(r" log $\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid("on")
p.tight_layout()

p.savefig(paper1_figures_path("M33_Sigma_profile_hi_co21_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(paper1_figures_path("M33_Sigma_profile_hi_co21_dr_{}pc.png".format(int(dr.value))))
p.close()

# p.show()

# Save the radial profiles in a FITS table
table = Table([rs, sd, sd_sigma, sd_hi, sd_sigma_hi],
              names=('Radius', "CO_Sigma", "CO_Sigma_std", "HI_Sigma",
                     "HI_Sigma_std"))
table.write(fourteenB_HI_data_path("tables/co21_hi_radialprofiles_{}pc.fits".format(int(dr.value)),
                                   no_check=True))

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
p.xlabel(r"Radius (kpc)")
p.legend(loc='best')
p.grid()
p.tight_layout()

p.savefig(paper1_figures_path("M33_Sigma_profile_gas_stars_corbelli_{}pc.pdf".format(int(dr.value))))
p.savefig(paper1_figures_path("M33_Sigma_profile_gas_stars_corbelli_{}pc.png".format(int(dr.value))))
p.close()


# p.show()

# Finally, let's calculate some clumping factors a la Leroy+13
rs_m, sd_m, sd_sigma_m = surfdens_radial_profile(gal, cube=cube, mom0=mom0,
                                                 max_rad=6 * u.kpc, dr=dr,
                                                 weight_type='mass',
                                                 mass_conversion=co21_mass_conversion)
# Correct for beam efficiency
sd_m /= beam_eff
sd_sigma_m /= beam_eff

rs_hi_m, sd_hi_m, sd_sigma_hi_m = \
    surfdens_radial_profile(gal, cube=hi_cube,
                            mom0=proj,
                            max_rad=6 * u.kpc, dr=dr,
                            weight_type='mass',
                            beam=average_beams(hi_cube.beams))

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
p.ylim([0.25, 1.9])
p.xlim([-2.1, 1.4])
p.legend(loc='upper left')
p.savefig(paper1_figures_path("hi_co_area_weighted_vs_mass_weighted_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(paper1_figures_path("hi_co_area_weighted_vs_mass_weighted_dr_{}pc.png".format(int(dr.value))))
p.tight_layout()
p.close()
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
p.ylim([0.65, 1.0])
p.xlim([0.65, 0.9])
p.tight_layout()
# p.legend(loc='upper left')

p.savefig(paper1_figures_path("area_weighted_vs_mass_weighted_dr_{}pc.pdf".format(int(dr.value))))
p.savefig(paper1_figures_path("area_weighted_vs_mass_weighted_dr_{}pc.png".format(int(dr.value))))
p.close()

# p.show()


clump_co = sd_m / sd
clump_hi = sd_hi_m / sd_hi

p.ion()
