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
os.sys.path.insert(0,parentdir)

from HI_veldisp_profile import radial_profile

# Load data in
direc = "/media/eric/Data_3/M33/co21/"
cube = SpectralCube.read(os.path.join(direc, "m33.co21_iram.fits"))
# cube = cube.with_mask(cube > 0.1 * u.K)
cube = cube.with_mask(fits.getdata(os.path.join(direc, "m33.co21_new_assign_clfind.fits")).astype(bool))

lwidth_co = cube.linewidth_sigma()

direc_hi = "/home/eric/MyRAID/M33/14B-088/HI/full_imaging/"
lwidth_hi = fits.open(os.path.join(direc_hi,
                                   "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.lwidth.fits"))[0]
lwidth_hi = Projection(lwidth_hi.data, wcs=WCS(lwidth_hi.header), unit=u.m / u.s)

g = Galaxy("M33")

rs, sd, sd_sigma = radial_profile(g, lwidth_co, max_rad=6 * u.kpc)

sd = sd.to(u.km / u.s)
sd_sigma = sd_sigma.to(u.km / u.s)

p.errorbar(rs.value, sd.value,
           yerr=sd_sigma.value, fmt="-", color="b",
           drawstyle='steps-mid')
p.xlabel("R (kpc)")
p.ylabel("CO Velocity Dispersion (km/s)")
p.grid()
p.draw()

raw_input("Next plot?")

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