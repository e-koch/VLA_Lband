
'''
Calculate luminosity of the one OH maser detection.
'''

from spectral_cube import SpectralCube
from astropy.io import fits
import pyregion
from os.path import join as osjoin
from astropy.wcs.utils import proj_plane_pixel_area
import astropy.units as u
import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt

from constants import distance


path = "/mnt/bigdata/ekoch/VLA/14B-088/Lines/OH/OH1665/imaging_1point5km_s"

# Open the OH 1665 narrow line width cube

cube = SpectralCube.read(osjoin(path, "OH1665_14B-088_uniform.image.pbcor.fits"))

reg = pyregion.open("/mnt/bigdata/ekoch/VLA/14B-088/Lines/OH/oh1665_maser_uniform.reg")

reg_cube = cube.subcube_from_ds9region(reg)
reg_cube_hz = reg_cube.with_spectral_unit(u.Hz, velocity_convention='radio', rest_value=1665 * u.MHz)
mom0 = reg_cube.moment0()
mom0_hz = reg_cube_hz.moment0()

# Integrated intensity over region

total_intensity = np.nansum(mom0).quantity
total_intensity_hz = np.nansum(mom0_hz).quantity

# Remove the per beam in the unit
num_pix = np.sum(np.isfinite(mom0))

pix_size = proj_plane_pixel_area(mom0.wcs) * u.Unit(mom0.header['CUNIT2'])**2

region_area = pix_size * num_pix

total_flux = total_intensity * cube.beam.sr.to(u.deg**2) * u.beam / region_area
total_flux_hz = total_intensity_hz * cube.beam.sr.to(u.deg**2) * u.beam / region_area

# Convert to SI units

total_flux_si = total_flux_hz.to(u.W / u.m**2)

total_luminosity = \
    4 * np.pi * (total_flux_si * distance.to(u.m)**2).to(u.solLum)
print("Maser luminosity of : {}".format(total_luminosity))
# Maser luminosity of : 0.000844790924261 solLum

# Per-channel limit to compare with Fix & Mutel 1985
# They list things in terms of Jy kpc^2. So this is a per channel limit?
# Take the peak intensity in the cube
peak_intensity = reg_cube.max()
peak_flux = peak_intensity * u.beam
peak_lum = 4 * np.pi * peak_flux * distance**2

print("Peak luminosity: {}".format(peak_lum))
# Peak luminosity: 163810.309843 Jy kpc2
# Fix & Mutel has a detection threshold of 1.3e4 Jy kpc^2, and this exceeds it
# Variability or dropping 4 pi?

# Comparison to some numbers from Walsh+16 for OH masers in W43.
# The brightest peak is
w43_peakflux = 24.395 * u.Jy
# Kinematic from Heyer+09
w43_dist = 5.7 * u.kpc

w43_chan_width = 0.73 * u.km / u.s
chan_width = np.diff(cube.spectral_axis[:2])[0]
chan_ratio = (w43_chan_width / chan_width.to(u.km / u.s))

w43_luminosity = w43_peakflux * w43_dist**2

lum_ratio = (peak_flux * distance**2 / w43_luminosity) * chan_ratio

print("Ratio b/w M33 peak luminosity and brightest in W43: {}"
      .format(lum_ratio))
# Ratio b/w M33 peak luminosity and brightest in W43: 8.0040986937
# So this maser is about ~10 times the brightest in W43

# Quick-and-dirty search for more emission.
# load in the flux coverage. This is a whole cube, but the coverage won't
# change drastically across different channels
flux = \
    fits.open(osjoin(path, "OH1665_14B-088_uniform.flux.pbcoverage.fits"))[0]
flux_mask = flux.data[0] > 0.7

cube = cube.with_mask(flux_mask)

# Calculate a threshold
mad_std = cube.mad_std()
# This doesn't seem to do a great job. It's quite high relative to the one
# detection. Look for things above this level as a first check.


# MAD std cut
masked_cube = cube.with_mask(cube > mad_std)

masked_mom0 = masked_cube.moment0()

# Calculate S/N ratio
noise_mom0 = \
    (mad_std * masked_cube.mask.include().sum(axis=0) / flux.data[0]) * \
    chan_width

# Slice out along flux cut
slice0 = nd.find_objects(np.isfinite(masked_cube[0]))[0]

snr_mom0 = masked_mom0[slice0].value / noise_mom0[slice0].value

snr_mask = snr_mom0 >= 3.

labels, num = nd.label(snr_mask, np.ones((3, 3)))

num_pix_labels = nd.sum(snr_mask, labels, range(1, num + 1))

# The beam is ~25 pixels. Allow regions half the size or so
for i in range(num):
    if num_pix_labels[i] < 13:
        snr_mask[labels == i + 1] = 0


# plt.imshow(snr_mom0, origin='lower')
# plt.colorbar()
# plt.contour(snr_mask, colors='r')

# And there is only one detection.

# What's our luminosity threshold?

peak_sigma = mad_std * u.beam
lum_sigma = peak_sigma * distance**2
print("Luminosity at 1-sigma: {}".format(lum_sigma))
# Luminosity at 1-sigma: 1232.20014071 Jy kpc2
# Pretty much the same as Fix & Mutel
# Not sure that makes sense...
