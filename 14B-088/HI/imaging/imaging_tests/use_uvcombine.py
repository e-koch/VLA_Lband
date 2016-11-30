
import numpy as np
import os
from radio_beam import Beam
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from image_tools.radialprofile import azimuthalAverage

from uvcombine.uvcombine import feather_kernel, fftmerge, feather_compare


data_path = os.path.expanduser("~/MyRAID/M33/VLA/14B-088/HI/feather_testing")


# high_hdu = fits.open(os.path.join(data_path, "M33_14B-088_HI.clean_channel_820.image_old.fits"))[0]
# high_hdu = fits.open(os.path.join(data_path, "14B-088_HI_LSRK.ms.contsub_channel_1514.clean.image.fits"))[0]
high_hdu = fits.open(os.path.join(data_path, "M33_14B-088_HI.clean_channel_820.image.fits"))[0]
# high_hdu_casa_feather = \
#     fits.open(os.path.join(data_path, "M33_14B-088_HI.clean.image.feathered_channel_820.image.fits"))[0]
low_hdu = fits.open(os.path.join(data_path, "M33_14B-088_HI_model_channel_844.image.fits"))[0]
mask = fits.open(os.path.join(data_path, "M33_14B-088_HI_mask_channel_844.fits"))[0]

high_beam = Beam.from_fits_header(high_hdu.header)
low_beam = Beam.from_fits_header(low_hdu.header)

y_size, x_size = high_hdu.shape

pixscale = (3 * u.arcsec).to(u.deg)

# Try removing all NaNed regions
# slicer = (slice(845, 2133), slice(725, 1615))
slicer = (slice(None), slice(None))

# new_shape = (2133 - 845, 1615 - 725)
new_shape = (2560, 2560)

# Ensure Arecibo doesn't have negative values or power outside of the VLA map
high_data = high_hdu.data.copy()[slicer]
low_data = low_hdu.data.copy()[::-1, ::-1][slicer]

low_data[low_data < 0.0] = 0.0
low_data[np.isnan(high_data)] = np.NaN

# kfft, ikfft = feather_kernel(2560, 2560, low_beam.major, pixscale)
kfft, ikfft = feather_kernel(new_shape[0], new_shape[1], low_beam.major, pixscale)

kernel_low = low_beam.as_kernel(pixscale, x_size=x_size, y_size=y_size)

# kfft = np.fft.fft2(kernel_low)
# ikfft = 1 - kfft

# Convert to Jy/pixel in both.
sr_pix = (pixscale ** 2).to(u.sr) / u.pix

low_data *= sr_pix / low_beam.sr
high_data *= sr_pix / high_beam.sr


fft_sum, merged = fftmerge(kfft, ikfft, high_data,
                           low_data,
                           highpassfilterSD=False,
                           replace_hires=False,
                           deconvSD=True, min_beam_fraction=0.1)

merged = merged.real

plt.imshow(merged, origin='lower')

# Compare power spectra.

high_fft = np.fft.fft2(np.nan_to_num(high_data))
low_fft = np.fft.fft2(np.nan_to_num(low_data))


# merged[np.isnan(high_hdu.data)] = np.NaN
merged_fft = np.fft.fft2(np.nan_to_num(merged))

# And the CASA feather version
# casa_fft = np.fft.fft2(np.nan_to_num(high_hdu_casa_feather.data))

rad, azavg_sum = azimuthalAverage(np.abs(np.fft.fftshift(merged_fft)),
                                  returnradii=True)

azavg_high = azimuthalAverage(np.abs(np.fft.fftshift(high_fft)))
azavg_low = azimuthalAverage(np.abs(np.fft.fftshift(low_fft)))

# azavg_casa = azimuthalAverage(np.abs(np.fft.fftshift(casa_fft)))

rad_pix = 2560. / rad
xaxis = pixscale.to(u.arcsec).value * rad_pix

ax1 = plt.subplot(1, 1, 1)
ax1.loglog(xaxis, azavg_high, color='b', linewidth=3,
           alpha=0.5,
           linestyle='-.',
           label="VLA")
ax1.loglog(xaxis, azavg_low, color='r', linewidth=2,
           alpha=0.5,
           linestyle='--',
           label="Arecibo")
ax1.loglog(xaxis, azavg_sum, color='g', linewidth=2,
           alpha=0.5,
           linestyle='--',
           label="Merged")
# ax1.loglog(xaxis, azavg_casa, color='k', linewidth=2,
#            alpha=0.5,
#            linestyle='--',
#            label="CASA Merged")

ax1.set_ylabel("Power spectrum $|FT|$")

ax1.grid(True)
ax1.legend(loc='lower right', frameon=True)
