
from astropy.io import fits
import numpy as np
from radio_beam import Beam
import astropy.units as u
from decimal import Decimal

'''
Total mass of HI in M33 from the zeroth moment.
'''


mom0 = fits.open("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits")[0]

beam = Beam.from_fits_header(mom0.header)

# In Mpc
distance = 0.84

# Several steps:.
# mom0 is in Jy/bm m/s, convert to Jy km/s
int_line = np.nansum(mom0.data) * 0.001 * \
    ((mom0.header["CDELT2"]**2) / (beam.sr.to(u.deg**2).value))
# Correction factor in the VLA data
int_line /= 1.45
# Conversion to mass
mass = 2.36e5 * (distance**2) * int_line
print("Total HI mass: {:.2E}".format(Decimal(mass)))
