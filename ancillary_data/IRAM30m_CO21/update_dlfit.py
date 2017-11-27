
'''
Update the CO21 and HI intensities in the dlfit.fits file, which contains fits
to the Draine models for M33. The old version used an older and incomplete IRAM
CO(2-1) map, and the HI are from the 14B-088 integrated intensity.
'''

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.utils.console import ProgressBar
from spectral_cube import Projection
from radio_beam import Beam
from os.path import join as osjoin

from paths import iram_co21_data_path, data_path, fourteenB_wGBT_HI_file_dict
from constants import hi_freq, beam_eff_30m_druard

tab = Table.read(osjoin(data_path, "dlfit.fits"))

co21_mom0 = Projection.from_hdu(fits.open(iram_co21_data_path("m33.co21_iram.mom0.fits"))[0])
co21_rms = Projection.from_hdu(fits.open(iram_co21_data_path("m33.rms.masked.fits"))[0])

hi_mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0])


# Convolve to the lowest resolution used for the fits: 500 um Herschel
# Actually @low-sky used a 60'' beam for the Draine model fitting
beam = Beam(60. * u.arcsec)

smooth_co21 = co21_mom0.convolve_to(beam)
# Remove regions outside of the original map extent
smooth_co21[np.isnan(co21_rms)] = np.NaN

# And the HI
smooth_hi = hi_mom0.convolve_to(beam)

# Now convert HI from Jy m / s to K km / s
smooth_hi = (smooth_hi.value * beam.jtok(hi_freq) / 1000.) * u.km / u.s

# Convert to K km/s
smooth_co21 = smooth_co21.to(u.K * u.km / u.s)

# Now match the points in the table to the map
dec_map, ra_map = smooth_co21.spatial_coordinate_map

spatial_coords = SkyCoord(ra=ra_map,
                          dec=dec_map,
                          frame='icrs')

dec_map, ra_map = hi_mom0.spatial_coordinate_map

spatial_coords_hi = SkyCoord(ra=ra_map,
                             dec=dec_map,
                             frame='icrs')

ra_map[np.isnan(co21_rms)] = np.NaN
dec_map[np.isnan(co21_rms)] = np.NaN

ra_limits = [np.nanmin(ra_map).value, np.nanmax(ra_map).value]
dec_limits = [np.nanmin(dec_map).value, np.nanmax(dec_map).value]


new_co21 = np.zeros((len(tab)))
new_hi = np.zeros((len(tab)))
for i, (ra, dec) in enumerate(ProgressBar(zip(tab['RA'], tab['DEC']))):

    posn = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')

    if ra < ra_limits[0] or ra > ra_limits[1] or dec < dec_limits[0] or dec > dec_limits[1]:
        new_co21[i] = np.NaN
    else:

        # min_posn = spatial_coords.separation(posn).argmin()
        mask = spatial_coords.separation(posn) < 60. * u.arcsec

        # twod_posn = np.unravel_index(min_posn, co21_mom0.shape)

        # new_co21[i] = smooth_co21[twod_posn].value
        new_co21[i] = np.nanmean(smooth_co21[mask].value)

    mask_hi = spatial_coords_hi.separation(posn) < 60. * u.arcsec

    new_hi[i] = np.nanmean(smooth_hi[mask_hi].value)


# Correct CO for beam efficiency
beam_eff = beam_eff_30m_druard
tab['CO21'] = Column(new_co21 / beam_eff)
tab['HI'] = Column(new_hi)

tab.write(osjoin(data_path, "updated_dlfit.fits"), format='fits',
          overwrite=True)
