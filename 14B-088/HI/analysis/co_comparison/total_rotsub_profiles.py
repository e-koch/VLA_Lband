
from astropy.utils.console import ProgressBar
import astropy.units as u
from spectral_cube import SpectralCube
import os
import numpy as np
import matplotlib.pyplot as p
from astropy.modeling import models, fitting

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

raw_input("Next plot?")

# Total CO mass. Using 6.7 Msol / pc^2 / K km s^-1\
pixscale = 0.84e6 * (np.pi / 180.) * np.abs(co_cube.header["CDELT2"]) * u.pc
chan_width = np.abs(co_cube.spectral_axis[1] - co_cube.spectral_axis[0]).to(u.km / u.s)
co21_mass_conversion = 6.7 * (u.Msun / u.pc ** 2) / (u.K * u.km / u.s)
beam_eff = 0.75  # Beam efficiency of IRAM @ 235 GHz
# Where total_spectrum_co is in K
total_co_mass = \
    total_spectrum_co[total_spectrum_co > 0 * u.K].sum() * chan_width * \
    pixscale ** 2 * co21_mass_conversion / beam_eff
print("Total H2 Mass from CO is {} Msol".format(total_co_mass))


# We want to find the widths of these profiles.
# CO is close enough to a gaussian
# The HI has wings. Are they real wings? Maybe.. but at least some portion
# is due to double-peaked spectra. We're going to model w/ two gaussians
# with the same mean.

# HI model
g_HI_init = models.Gaussian1D(amplitude=1., mean=0., stddev=5.) +  \
    models.Gaussian1D(amplitude=0.25, mean=0., stddev=20.)

# Force to the same mean
def tie_mean(mod):
    return mod.mean_0
g_HI_init.mean_1.tied = tie_mean

fit_g = fitting.LevMarLSQFitter()

vels = hi_cube.spectral_axis.to(u.km / u.s).value
norm_intens = (total_spectrum_hi / total_spectrum_hi.max()).value
g_HI = fit_g(g_HI_init, vels, norm_intens)

# The covariance matrix is hidden away... tricksy
cov = fit_g.fit_info['param_cov']
parnames = [n for n in g_HI.param_names if n not in ['mean_1']]
parvals = [v for (n, v) in zip(g_HI.param_names, g_HI.parameters)
           if n in parnames]
print("HI Fit")
for i, (name, value) in enumerate(zip(parnames, parvals)):
    print('{}: {} +/- {}'.format(name, value, np.sqrt(cov[i][i])))

# Note that the statistical errors on the mean are too small.
# Due to my rolling approximation, the error is ~1 channel width.

p.plot(vels, norm_intens, 'b-', drawstyle='steps-mid')
p.plot(vels, g_HI(vels), 'k:', label="Total Fit")
p.plot(vels, g_HI["None_0"](vels), 'g--', label="Narrow Component")
p.plot(vels, g_HI["None_1"](vels), 'm-.', label="Wide Component")
p.xlabel("Velocity (km/s)")
p.ylabel("Total Normalized Intensity")
p.xlim([-100, 100])
p.legend()
p.ylim([-0.1, 1.1])
p.grid()
p.show()

raw_input("Next plot?")

# CO model

g_CO_init = models.Gaussian1D(amplitude=1., mean=0., stddev=9.)

fit_g_co = fitting.LevMarLSQFitter()

co_vels = co_cube.spectral_axis.to(u.km / u.s).value
norm_co_intens = np.roll((total_spectrum_co / total_spectrum_co.max()).value,
                         -1)
g_CO = fit_g_co(g_CO_init, co_vels, norm_co_intens)

cov = fit_g.fit_info['param_cov']
print("CO(2-1) Fit")
for i, (name, value) in enumerate(zip(g_CO.param_names, g_CO.parameters)):
    print('{}: {} +/- {}'.format(name, value, np.sqrt(cov[i][i])))

# Better sampling for plotting
more_vels = np.arange(co_vels.min(), co_vels.max(), 0.5)

p.plot(co_vels, norm_co_intens, 'b-', drawstyle='steps-mid')
p.plot(more_vels, g_CO(more_vels), 'k--', label="Total Fit")
p.xlabel("Velocity (km/s)")
p.ylabel("Total Normalized Intensity")
p.xlim([-100, 100])
p.ylim([-0.1, 1.1])
p.grid()
p.show()