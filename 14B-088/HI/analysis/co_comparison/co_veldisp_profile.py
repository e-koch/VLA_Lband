from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from astropy.wcs import WCS
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle

from galaxies import Galaxy
import matplotlib.pyplot as p
import os
from astropy.io import fits


# Import from above.
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.insert(0, parentdir)

from paths import (iram_co21_data_path, fourteenB_HI_data_path,
                   paper1_figures_path)
from constants import lwidth_name

from HI_veldisp_profile import radial_profile

# Load data in
cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
# cube = cube.with_mask(cube > 0.1 * u.K)
mask = fits.getdata(iram_co21_data_path("m33.co21_new_assign_clfind.fits"))
cube = cube.with_mask(mask.astype(bool))

lwidth_co = cube.linewidth_sigma()

lwidth_hi = fits.open(fourteenB_HI_data_path(lwidth_name))[0]
lwidth_hi = Projection(lwidth_hi.data, wcs=WCS(lwidth_hi.header),
                       unit=u.m / u.s)

g = Galaxy("M33")

dr = 250 * u.pc

rs, sd, sd_sigma = radial_profile(g, lwidth_co, max_rad=6 * u.kpc,
                                  dr=dr)

sd = sd.to(u.km / u.s)
sd_sigma = sd_sigma.to(u.km / u.s)

p.errorbar(rs.value, sd.value,
           yerr=sd_sigma.value, fmt="-", color="b",
           drawstyle='steps-mid')
p.xlabel("R (kpc)")
p.ylabel("CO Velocity Dispersion (km/s)")
p.grid()
p.draw()

p.savefig(paper1_figures_path("co_veldisp_profile_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("co_veldisp_profile_{}pc.png".format(dr.value)))
p.close()

# raw_input("Next plot?")

# Create the North and South portions.
rs_n, sd_n, sd_sigma_n = \
    radial_profile(g, lwidth_co, max_rad=6 * u.kpc,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
sd_n = sd_n.to(u.km / u.s)
sd_sigma_n = sd_sigma_n.to(u.km / u.s)
rs_s, sd_s, sd_sigma_s = \
    radial_profile(g, lwidth_co, max_rad=6 * u.kpc,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))
sd_s = sd_s.to(u.km / u.s)
sd_sigma_s = sd_sigma_s.to(u.km / u.s)

p.plot(rs.value, sd.value, "k-.",
       drawstyle='steps-mid', label="Total")
p.plot(rs_n.value, sd_n.value, "b-", label="North",
       drawstyle='steps-mid')
p.plot(rs_s.value, sd_s.value, "g--", label="South",
       drawstyle='steps-mid')
p.xlabel("R (kpc)")
p.ylabel("CO Velocity Dispersion (km/s)")
p.grid()
p.legend()
p.draw()

p.savefig(paper1_figures_path("co_veldisp_profile_n_s_{}pc.pdf".format(dr.value)))
p.savefig(paper1_figures_path("co_veldisp_profile_n_s_{}pc.png".format(dr.value)))
p.close()
