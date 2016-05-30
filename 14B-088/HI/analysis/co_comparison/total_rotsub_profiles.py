
from astropy.utils.console import ProgressBar
import astropy.units as u
from spectral_cube import SpectralCube
import os
import numpy as np
import matplotlib.pyplot as p

'''
Create profiles of HI and CO after rotation subtraction.
'''

co_data_path = "/media/eric/Data_3/M33/co21/"
hi_data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"


co_cube = SpectralCube.read(os.path.join(co_data_path,
                            "m33.co21_iram.rotsub.fits"))

hi_cube = SpectralCube.read(os.path.join(hi_data_path,
                            "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits"))

total_spectrum_hi = np.empty((hi_cube.shape[0]))
total_spectrum_co = np.empty((co_cube.shape[0]))

# Create the HI spectrum
for chan in ProgressBar(range(hi_cube.shape[0])):
    channel = hi_cube[chan]
    beam = channel.meta['beam']
    # Values are in Jy/bm. To get Jy, convert from bm to deg**2, then multiply
    # by the pixel area in deg^2
    total_spectrum_hi[chan] = \
        np.nansum(channel.value) * (1 / beam.sr.to(u.deg ** 2)) * \
        (channel.header["CDELT2"] * u.deg)**2
    # total_spectrum_hi[chan] = (total_spectrum_hi[chan] * u.Jy).to(u.K, equivalencies=beam.jtok_equiv(1.4 * u.GHz)).value
total_spectrum_hi
for chan in ProgressBar(range(co_cube.shape[0])):
    channel = co_cube[chan]
    beam = channel.meta['beam']
    # Units are in K, but applying the pixel area factor early. Jy conversion below
    total_spectrum_co[chan] = \
        np.nansum(channel.value) # * (1 / beam.sr.to(u.deg ** 2)) * \
        # (channel.header["CDELT2"] * u.deg)**2
    # CO cube is in K. Convert to Jy
    # total_spectrum_co[chan] = \
    #     (total_spectrum_co[chan] * u.K).to(u.Jy, equivalencies=beam.jtok_equiv(230.538 * u.GHz)).value

# VLA needs to be corrected by a 1.45 factor. See arecibo_match_props.py
total_spectrum_hi = total_spectrum_hi * u.Jy / 1.45
total_spectrum_co = total_spectrum_co * u.K


p.plot(hi_cube.spectral_axis.to(u.km / u.s).value, total_spectrum_hi.value,
       'b-', drawstyle='steps-mid', label="HI")
p.xlabel("Velocity (km/s)")
p.ylabel("Total Intensity (Jy)")
p.xlim([-100, 100])
p.show()

raw_input("Next plot?")
p.clf()

p.plot(hi_cube.spectral_axis.to(u.km / u.s).value,
       (total_spectrum_hi / total_spectrum_hi.max()).value,
       'b-', drawstyle='steps-mid', label="HI")
# There's a 1 channel offset from my rotation subtraction in the cube
p.plot(co_cube.spectral_axis.to(u.km / u.s).value,
       np.roll((total_spectrum_co / total_spectrum_co.max()).value, -1),
       'g--', drawstyle='steps-mid', label="CO(2-1)")
p.xlabel("Velocity (km/s)")
p.ylabel("Normalized Total Intensity")
p.ylim([-0.02, 1.1])
p.xlim([-100, 100])
p.legend()
p.show()

# Total CO mass. Using 6.7 Msol / pc^2 / K km s^-1\
pixscale = 0.84e6 * (np.pi / 180.) * np.abs(co_cube.header["CDELT2"]) * u.pc
chan_width = np.abs(co_cube.spectral_axis[1] - co_cube.spectral_axis[0]).to(u.km / u.s)
co21_mass_conversion = 6.7 * (u.Msun / u.pc ** 2) / (u.K * u.km / u.s)
beam_eff = 0.75  # Beam efficiency of IRAM @ 235 GHz
# Where total_spectrum_co is in K
total_co_mass = \
    total_spectrum_co[total_spectrum_co > 0 * u.K].sum() * chan_width * \
    pixscale ** 2 * co21_mass_conversion / beam_eff
