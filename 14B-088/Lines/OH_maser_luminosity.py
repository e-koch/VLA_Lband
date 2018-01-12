
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


path = "/mnt/bigdata/ekoch/VLA/14B-088/Lines/OH/OH1665/detection_imaging_1point5km_s"

# Open the OH 1665 narrow line width cube

cube = SpectralCube.read(osjoin(path, "OH1665_14B-088_uniform.image.pbcor.fits"))

# Convolve to a common beam
com_beam = cube.beams.common_beam(epsilon=7e-4)
cube = cube.convolve_to(com_beam)

reg = pyregion.open("/mnt/bigdata/ekoch/VLA/14B-088/Lines/OH/oh1665_maser_uniform.reg")

reg_cube = cube.subcube_from_ds9region(reg)
reg_cube_hz = reg_cube.spectral_slab(-220 * u.km / u.s, -200 * u.km / u.s)\
    .with_spectral_unit(u.Hz, velocity_convention='radio',
                        rest_value=1665 * u.MHz)
mom0 = reg_cube.spectral_slab(-220 * u.km / u.s, -200 * u.km / u.s).moment0()
mom0_hz = reg_cube_hz.moment0()

# Integrated intensity over region

num_pix = np.sum(np.isfinite(mom0)) * u.pix
pix_size = proj_plane_pixel_area(mom0.wcs) * u.Unit(mom0.header['CUNIT2'])**2

pix_per_beam = (com_beam.sr.to(u.deg**2) / pix_size) * (u.pix / u.beam)

total_flux = np.nansum(mom0).quantity * (num_pix / pix_per_beam)
total_flux_hz = np.nansum(mom0_hz).quantity * (num_pix / pix_per_beam)

# Make a total spectrum
sum_spec = reg_cube.sum(axis=(1, 2)) * (num_pix / pix_per_beam)

# Convert to SI units

total_flux_si = total_flux_hz.to(u.W / u.m**2)

total_luminosity = \
    4 * np.pi * (total_flux_si * distance.to(u.m)**2).to(u.solLum)
print("Maser luminosity of : {}".format(total_luminosity))
# Maser luminosity of : 0.00618556389833 solLum
print(total_luminosity.to(u.erg / u.s))

# Per-channel limit to compare with Fix & Mutel 1985
# They list things in terms of Jy kpc^2. So this is a per channel limit?
# Take the peak intensity in the summed spectrum
peak_flux = sum_spec.max()
peak_lum = 4 * np.pi * peak_flux * distance**2

print("Peak luminosity: {}".format(peak_lum))
# Peak luminosity: 8219629.01471 Jy kpc2

# Fix & Mutel has a detection threshold of 1.3e4 Jy kpc^2 @ 720 kpc distance
fix_mutel_5sig_lim = 1.3e4 * u.Jy * u.kpc**2 * (840 / 720.)**2
# @ 840 kpc, 1.77e4 Jy kpc^2
# This is NOT an isotropic value. Make it isotropic to compare with our peak
fix_mutel_5sig_lim_iso = fix_mutel_5sig_lim * 4 * np.pi

# And finally correct for the difference in channel size. They have 5.2 km/s
fix_mutel_5sig_matched = fix_mutel_5sig_lim_iso * np.sqrt(5.2 / 1.5)

# So the S/N in the Fix & Mutel data would be:
print("Fix & Mutel S/N would be: {}".format(5 * (peak_lum / fix_mutel_5sig_matched)))
# Fix & Mutel S/N would be: 99.2703111135
# It SHOULD have been detectable!

# Comparison to some numbers from Walsh+16 for OH masers in W43.
# The brightest peak is
w43_peakflux = 24.395 * u.Jy
# Kinematic from Heyer+09
w43_dist = 5.5 * u.kpc

w43_chan_width = 0.73 * u.km / u.s
chan_width = np.diff(cube.spectral_axis[:2])[0]
chan_ratio = (w43_chan_width / chan_width.to(u.km / u.s))

w43_luminosity = 4 * np.pi * w43_peakflux * w43_dist**2

lum_ratio = (peak_lum / w43_luminosity) * np.sqrt(chan_ratio)

print("Ratio b/w M33 peak luminosity and brightest in W43: {}"
      .format(lum_ratio))
# Ratio b/w M33 peak luminosity and brightest in W43: 618.346526462
# So this maser is about ~600 times the brightest in W43

# Quick-and-dirty search for more emission.
# load in the flux coverage. This is a whole cube, but the coverage won't
# change drastically across different channels

# I originally ran this on the whole map, and didn't find any sources. The
# following just confirms the one detection is a good detection

flux = fits.open(osjoin(path, "OH1665_14B-088_uniform.pb.fits"))[0]
flux_mask = flux.data[0] > 0.7

cube = cube.with_mask(flux_mask)

# Calculate a threshold
mad_std = cube.mad_std()
# Using a large circular region in casaview, I get a std of ~1.3e-3 Jy/bm
# The mad_std gives 1.86 mJy/bm. Accounting for pb variations, this seems
# reasonable

# MAD std cut
masked_cube = cube.with_mask(cube > mad_std).spectral_slab(-220 * u.km / u.s, -200 * u.km / u.s)

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
    if num_pix_labels[i] < 10:
        snr_mask[labels == i + 1] = 0


# plt.imshow(snr_mom0, origin='lower')
# plt.colorbar()
# plt.contour(snr_mask, colors='r')

# And there is only one detection.

# What's our luminosity threshold?

peak_sigma = mad_std * (num_pix / pix_per_beam)
lum_sigma = 4 * np.pi * peak_sigma * distance**2
print("Luminosity at 1-sigma: {}".format(lum_sigma))
# Luminosity at 1-sigma: 40202.6195055 Jy kpc2