
'''
Make a table of azimuthal avgs in radial bins for several
tracers.
'''

from astropy.io import fits
from spectral_cube import Projection
import numpy as np
import astropy.units as u
from os.path import join as osjoin
import os
from astropy.table import Table, Column

from cube_analysis.profiles import radial_profile

from paths import (data_path, fourteenB_wGBT_HI_file_dict,
                   fourteenB_HI_data_wGBT_path,
                   iram_co21_14B088_data_path, allfigs_path)
# from paths import data_path
from constants import (co21_mass_conversion, hi_mass_conversion, hi_freq,
                       beam_eff_30m_druard)
from plotting_styles import default_figure, onecolumn_figure, twocolumn_figure
from galaxy_params import gal_feath as gal


dust_col_hdu = fits.open(osjoin(data_path, "m33_dust.surface.density_FB.beta=1.8_gauss41.0_regrid_bksub.fits"))[0]
dust_temp_hdu = fits.open(osjoin(data_path, "m33_dust.temperature_FB.beta=1.8_gauss41.0_regrid_bksub.fits"))[0]

hi_mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0']))
co_mom0 = Projection.from_hdu(fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits")))

# Convert to K km/s in HI
hi_mom0_Kkms = (hi_mom0 * hi_mom0.beam.jtok(hi_freq) * u.beam / u.Jy).to(u.K * u.km / u.s) * np.cos(gal.inclination)

hi_surfdens = hi_mass_conversion * hi_mom0_Kkms

co_mom0 = co_mom0.to(u.K * u.km / u.s) * np.cos(gal.inclination) / beam_eff_30m_druard
co_surfdens = co_mom0 * co21_mass_conversion

co_mask = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_source_mask.fits"))[0]

# good_pts = np.where(np.isfinite(mom0_reproj))
good_co_mask = co_mask.data.sum(0) >= 2

co_surfdens = co_surfdens.with_mask(good_co_mask)

# Radius
radius = gal.radius(header=hi_mom0.header).to(u.kpc)

dust_col_hdr = dust_col_hdu.header.copy()
dust_col_hdr['BUNIT'] = "solMass / pc2"
dust_surfdens = Projection.from_hdu(fits.PrimaryHDU(dust_col_hdu.data[0],
                                                    dust_col_hdr))
dust_surfdens_rep = dust_surfdens.reproject(hi_mom0.header)

dust_temp_hdr = dust_temp_hdu.header.copy()
dust_temp_hdr['BUNIT'] = "K"
dust_temp = Projection.from_hdu(fits.PrimaryHDU(dust_temp_hdu.data[0],
                                                dust_temp_hdr))
dust_temp_rep = dust_temp.reproject(hi_mom0.header)

stellar_surfdens = Projection.from_hdu(fits.open(osjoin(data_path, "m33.stellarmass.fits")))
stellar_surfdens = stellar_surfdens * u.solMass / u.kpc**2
stellar_surfdens_rep = stellar_surfdens.reproject(hi_mom0.header).to(u.solMass / u.pc**2)

# Make avgs of surface densities and of observed quantities (where available)

dr = 500 * u.pc

rs, hi_intint_avg, hi_intint_sigma = \
    radial_profile(gal, hi_mom0_Kkms, max_rad=6 * u.kpc,
                   dr=dr)

# Should all have the same bins
radius = Column(rs, name='Radius')
hi_intint_avg = Column(hi_intint_avg, name='hi_intint')
hi_intint_sigma = Column(hi_intint_sigma, name='hi_intint_sigma')

rs, hi_sd_avg, hi_sd_sigma = \
    radial_profile(gal, hi_surfdens, max_rad=6 * u.kpc,
                   dr=dr)
hi_sd_avg = Column(hi_sd_avg, name='hi_sd')
hi_sd_sigma = Column(hi_sd_sigma, name='hi_sd_sigma')

# Dust

rs, dust_sd_avg, dust_sd_sigma = \
    radial_profile(gal, dust_surfdens_rep, max_rad=6 * u.kpc,
                   dr=dr)
dust_sd_avg = Column(dust_sd_avg, name='dust_sd')
dust_sd_sigma = Column(dust_sd_sigma, name='dust_sd_sigma')

# Stellar mass

rs, stellar_sd_avg, stellar_sd_sigma = \
    radial_profile(gal, stellar_surfdens_rep, max_rad=6 * u.kpc,
                   dr=dr)
stellar_sd_avg = Column(stellar_sd_avg, name='stellar_sd')
stellar_sd_sigma = Column(stellar_sd_sigma, name='stellar_sd_sigma')

# CO

rs, co_intint_avg, co_intint_sigma = \
    radial_profile(gal, co_mom0, max_rad=6 * u.kpc,
                   dr=dr)

# Should all have the same bins
co_intint_avg = Column(co_intint_avg, name='co_intint')
co_intint_sigma = Column(co_intint_sigma, name='co_intint_sigma')

rs, co_sd_avg, co_sd_sigma = \
    radial_profile(gal, co_surfdens, max_rad=6 * u.kpc,
                   dr=dr)
co_sd_avg = Column(co_sd_avg, name='co_sd')
co_sd_sigma = Column(co_sd_sigma, name='co_sd_sigma')

tab = Table([radius, hi_intint_avg, hi_intint_sigma, hi_sd_avg, hi_sd_sigma,
             co_intint_avg, co_intint_sigma, co_sd_avg, co_sd_sigma,
             dust_sd_avg, dust_sd_sigma, stellar_sd_avg, stellar_sd_sigma])
tab.write(fourteenB_HI_data_wGBT_path("tables/sd_radial_500pc_w_stellar_dust.fits",
                                      no_check=True),
          overwrite=True)
