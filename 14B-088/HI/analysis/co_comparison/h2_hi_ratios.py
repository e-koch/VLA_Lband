
'''
Create tables/maps of the H2/HI ratio and related quantities.
'''

import numpy as np
from astropy import units as u
import astropy.io.fits as fits
from astropy.table import Table
from spectral_cube import Projection
import os

from paths import (iram_co21_14B088_data_path, fourteenB_HI_data_path,
                   fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   fourteenB_HI_data_wGBT_path)
from constants import (co21_mass_conversion, hi_mass_conversion, hi_freq)
from galaxy_params import gal_feath as gal

# Now load in the zeroth moments.

hi_mom0 = Projection.from_hdu(fits.open(fourteenB_HI_file_dict['Moment0'])[0])
# Convert to K km/s
hi_mom0_data = hi_mom0.value * hi_mom0.beam.jtok(hi_freq).value / 1000.
hi_mom0_data = hi_mom0_data * u.K * u.km / u.s

mom0_reproj = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits"))[0]
mom0_reproj = (mom0_reproj.data / 1000.) * u.K * u.km / u.s

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

hi_coldens_cut = hi_coldens[~nans]
co_coldens_cut = co_coldens[~nans]
radii_pts_cut = radii_pts[~nans]
pang_pts_cut = pang_pts[~nans]
skycoord_pts_cut = skycoord_pts[~nans]
ypts_cut = good_pts[0][~nans]
xpts_cut = good_pts[1][~nans]

gas_ratio_pix = co_coldens_cut / hi_coldens_cut
total_sd_pix = co_coldens_cut + hi_coldens_cut

if not os.path.exists(fourteenB_HI_data_path("tables", no_check=True)):
    os.mkdir(fourteenB_HI_data_path("tables", no_check=True))

# Save the lists of points in a table
tab = Table([skycoord_pts_cut.ra, skycoord_pts_cut.dec, radii_pts_cut,
             pang_pts_cut,
             hi_coldens_cut, co_coldens_cut, total_sd_pix, gas_ratio_pix,
             ypts_cut, xpts_cut],
            names=["RA", "Dec", "Radius", "PA", "Sigma_HI", "Sigma_H2",
                   "Sigma_Total", "Ratio", "xpix", "ypix"])
tab.write(fourteenB_HI_data_path("tables/column_densities_perpix.fits",
                                 no_check=True),
          overwrite=True)

# Make another table with the feathered HI data.
hi_mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0])
# Convert to K km/s
hi_mom0_data = hi_mom0.value * hi_mom0.beam.jtok(hi_freq).value / 1000.
hi_mom0_data = hi_mom0_data * u.K * u.km / u.s
hi_coldens = hi_mom0_data[good_pts] * hi_mass_conversion * inc

nans = np.logical_or(np.isnan(hi_coldens), np.isnan(co_coldens))

hi_coldens_cut = hi_coldens[~nans]
co_coldens_cut = co_coldens[~nans]
radii_pts_cut = radii_pts[~nans]
pang_pts_cut = pang_pts[~nans]
skycoord_pts_cut = skycoord_pts[~nans]
ypts_cut = good_pts[0][~nans]
xpts_cut = good_pts[1][~nans]

gas_ratio_pix = co_coldens_cut / hi_coldens_cut
total_sd_pix = co_coldens_cut + hi_coldens_cut

if not os.path.exists(fourteenB_HI_data_wGBT_path("tables", no_check=True)):
    os.mkdir(fourteenB_HI_data_wGBT_path("tables", no_check=True))

# Save the lists of points in a table
tab = Table([skycoord_pts_cut.ra, skycoord_pts_cut.dec, radii_pts_cut,
             pang_pts_cut,
             hi_coldens_cut, co_coldens_cut, total_sd_pix, gas_ratio_pix,
             ypts_cut, xpts_cut],
            names=["RA", "Dec", "Radius", "PA", "Sigma_HI", "Sigma_H2",
                   "Sigma_Total", "Ratio", "xpix", "ypix"])
tab.write(fourteenB_HI_data_wGBT_path("tables/column_densities_perpix.fits",
                                      no_check=True),
          overwrite=True)
