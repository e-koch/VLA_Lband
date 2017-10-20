
'''
Masses and uncertainties from the zeroth moments.
'''

import astropy.units as u
import numpy as np
from astropy.io import fits
from spectral_cube import Projection
from pandas import DataFrame

from paths import (fourteenB_wGBT_HI_file_dict, fourteenB_HI_file_dict,
                   gbt_HI_data_path, alltables_path,
                   fourteenB_HI_data_wGBT_path)
from constants import (hi_mass_conversion_Jy, distance, hi_freq,
                       hi_mass_conversion)

mom0_feath = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0']))
mom0 = Projection.from_hdu(fits.open(fourteenB_HI_file_dict['Moment0']))

gbt_name = gbt_HI_data_path("14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid_registered.mom0.fits")
mom0_gbt = Projection.from_hdu(fits.open(gbt_name)[0])

feath_mask = fits.open(fourteenB_wGBT_HI_file_dict['Source_Mask'])[0]
mask = fits.open(fourteenB_HI_file_dict['Source_Mask'])[0]

feath_spatial_mask = feath_mask.data.sum(0)
spatial_mask = mask.data.sum(0)

sigma_noise = 2.8 * u.K

beam = mom0_feath.beam
pix_area_sr = ((mom0_feath.header['CDELT2'] * u.deg)**2).to(u.sr)
beams_per_pix = pix_area_sr / beam.sr

sigma_noise_Jybm = sigma_noise.to(u.Jy, beam.jtok_equiv(hi_freq))

vla_mass = np.nansum(mom0).quantity.to(u.Jy * u.km / u.s) * \
    hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

vla_mass_stddev = (sigma_noise_Jybm * 0.2 * u.km / u.s * spatial_mask).sum() * \
    hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

feath_mass = np.nansum(mom0_feath).quantity.to(u.Jy * u.km / u.s) * \
    hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

feath_mass_stddev = \
    (sigma_noise_Jybm * 0.2 * u.km / u.s * feath_spatial_mask).sum() * \
    hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * beams_per_pix

gbt_beam = mom0_gbt.beam
gbt_pix_area_sr = ((mom0_gbt.header['CDELT2'] * u.deg)**2).to(u.sr)
gbt_beams_per_pix = gbt_pix_area_sr / gbt_beam.sr

gbt_sigma_noise = 0.06 * u.K
gbt_sigma_noise_Jybm = gbt_sigma_noise.to(u.Jy, gbt_beam.jtok_equiv(hi_freq))

gbt_mom0_sum = np.nansum(mom0_gbt).quantity / gbt_beam.jtok(hi_freq) * u.Jy

gbt_mass = gbt_mom0_sum * hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * gbt_beams_per_pix

# No masking for the GBT cube, so use all channels
nchans = 1178

gbt_mass_stddev = \
    (gbt_sigma_noise_Jybm * 0.2 * u.km / u.s * 1178 * mom0_gbt.size) * \
    hi_mass_conversion_Jy * \
    distance.to(u.Mpc)**2 * gbt_beams_per_pix


print("VLA HI Total Mass: {0}+\-{1}".format(vla_mass, vla_mass_stddev))
print("GBT HI Total Mass: {0}+\-{1}".format(gbt_mass, gbt_mass_stddev))
print("VLA + GBT HI Total Mass: {0}+\-{1}".format(feath_mass,
                                                  feath_mass_stddev))

df = DataFrame({"VLA Mass": [vla_mass.value, vla_mass_stddev.value],
                "GBT Mass": [gbt_mass.value, gbt_mass_stddev.value],
                "VLA+GBT Mass": [feath_mass.value, feath_mass_stddev]})
df.to_latex(alltables_path("hi_masses_from_mom0.tex"))
df.to_csv(fourteenB_HI_data_wGBT_path("tables/hi_masses_from_mom0.csv",
                                      no_check=True))
