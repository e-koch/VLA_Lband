
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.utils.console import ProgressBar
import os
import matplotlib.pyplot as p
import numpy as np

'''
Sum over all spectra above 3 sigma. Overplot with single-dish Arecibo.
And produce the summed spectra of the rotation subtracted cube.
'''

# You know what's great? Making an install script!
execfile(os.path.expanduser("~/Dropbox/code_development/BaSiCs/basics/utils.py"))
# Needed sig_clip.


# p.ioff()

data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"
arecibo_path = "/media/eric/Data_3//M33/Arecibo/14B-088_items/"

cube = \
    SpectralCube.read(os.path.join(data_path,
                                   "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))
                                   # "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits"))
arecibo = SpectralCube.read(os.path.join(arecibo_path, "M33_14B-088_HI_model.fits"))
arecibo = arecibo.spectral_slab(*cube.spectral_extrema)
arecibo = arecibo.subcube(ylo=cube.latitude_extrema[1],
                          yhi=cube.latitude_extrema[0],
                          xlo=cube.longitude_extrema[0],
                          xhi=cube.longitude_extrema[1])

# noise_sigma = sig_clip(cube[-1].value, nsig=10) * u.Jy

# cube = cube.with_mask(cube > 3 * noise_sigma)

restfreq = 1.414 * u.GHz

# Using sum on two axes of a cube is super slow (and wrong? All NaNs?).
# I'll just loop through each channel to do this
total_spectrum = np.empty((cube.shape[0]))
total_spectrum_arecibo = np.empty((cube.shape[0]))
for chan in ProgressBar(range(cube.shape[0])):
    channel = cube[chan]
    beam = channel.meta['beam']
    # Values are in Jy/bm. To get Jy, convert from bm to deg**2, then multiply
    # by the pixel area in deg^2
    total_spectrum[chan] = \
        np.nansum(channel.value) * (1 / beam.sr.to(u.deg ** 2)) * \
        (channel.header["CDELT2"] * u.deg)**2

    arec_channel = arecibo[chan]
    arec_beam = arec_channel.meta['beam']
    total_spectrum_arecibo[chan] = \
        np.nansum(arec_channel.value) * (1 / arec_beam.sr.to(u.deg ** 2)) * \
        (arec_channel.header["CDELT2"] * u.deg)**2

# VLA needs to be corrected by a 1.45 factor. See arecibo_match_props.py
total_spectrum = total_spectrum * u.Jy / 1.45
total_spectrum_arecibo = total_spectrum_arecibo * u.Jy

p.plot(cube.spectral_axis.to(u.km / u.s).value, total_spectrum.value, 'b-',
       label="VLA + Arecibo", drawstyle='steps-mid')
p.plot(arecibo.spectral_axis.to(u.km / u.s).value,
       total_spectrum_arecibo.value, 'g--', label='Arecibo',
       drawstyle='steps-mid')
p.ylabel("Total Intensity (Jy)")
p.xlabel("Velocity (km/s)")
p.legend(loc='upper left')

# Now we can get the total mass

# Disregard outside this velocity range
vmax = -100 * u.km / u.s
vmin = -300 * u.km / u.s

spec_axis = cube.spectral_axis.to(u.km / u.s)
good_chans = np.logical_and(spec_axis <= vmax, spec_axis >= vmin)

# Same channel width in both
chan_width = \
    np.abs(cube.spectral_axis[1] - cube.spectral_axis[0]).to(u.km / u.s)

total_line_flux = np.nansum(total_spectrum[good_chans]) * chan_width
total_line_flux_arecibo = \
    np.nansum(total_spectrum_arecibo[good_chans]) * chan_width

# Distance is 0.84 Mpc
mass_conversion = 2.36e5 * (0.84 ** 2)
# These estimates are low by a factor of ~0.85. Not entirely sure why, but
# the estimate from the moment 0 is correct (comparing to Putman+09)

# total_mass = total_line_flux * mass_conversion
# total_mass_arecibo = total_line_flux_arecibo * mass_conversion
# print("VLA HI Total Mass: {} Msun".format(total_mass.value))
# print("Arecibo HI Total Mass: {} Msun".format(total_mass_arecibo.value))
