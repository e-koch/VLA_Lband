
'''
Use feather plot to examine differences in the cleaned image to Arecibo.
'''

import os
from glob import glob
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from radio_beam import Beam

from astropy.utils.console import ProgressBar
import astropy.io.fits as fits
from image_tools.radialprofile import azimuthalAverage

from uvcombine.feather_plot import feather_plot

data_path = os.path.expanduser("~/MyRAID/M33/14B-088/HI/channel_testing")
# data_path = os.path.expanduser("/Volumes/Travel_Data/M33/testing")


def append_path(path):
    return os.path.join(data_path, path)


parameters = ["CASAVer", "Model", "Mask", "AllFields", "MScale", "Tclean"]

# Load in the mask
mask = \
    fits.open(append_path("M33_14B-088_HI_mask_channel_330.fits"))[0].data.squeeze()
mask = mask[::-1, ::-1] == 1

model = \
    fits.open(append_path("M33_14B-088_HI_model_channel_330.fits"))[0].data.squeeze()
low_hdu = fits.open(append_path("M33_14B-088_HI_model_channel_330.fits"))[0]
# For some reason, both axes need to be reversed
model = model[::-1, ::-1]
sd_beam = Beam.from_fits_header(append_path("M33_14B-088_HI_model_channel_330.fits"))
fft_sd = np.fft.fftshift(np.fft.fft2(model))

# Arecibo FWHM
lowresfwhm = (3.8 * u.arcmin).to(u.arcsec)
# Grid scaling used.
pixscale = 3 * u.arcsec
# Image shape
nax2, nax1 = model.shape

# High scale factor
highresscalefactor = 1.0
# Low scale factor
lowresscalefactor = 1.0

# x-axis unit
xaxisunit = 'lambda'

sd_kernel = sd_beam.as_kernel(pixscale, x_size=nax2, y_size=nax1)
kfft = np.abs(np.fft.fft2(sd_kernel))
ikfft = 1 - kfft

rad, azavg_kernel = azimuthalAverage(np.abs(kfft), returnradii=True)
rad, azavg_ikernel = azimuthalAverage(np.abs(ikfft), returnradii=True)

fft_lo = \
    np.fft.fftshift(
        np.fft.fft2(
            np.nan_to_num(model * lowresscalefactor)))

folders = glob(append_path("14B-088_HI_LSRK.ms.contsub_channel_1000.CASAVer*"))

# Values to extract out.
for folder in folders:

    # Slice out the first 3 since this is just the name
    split_name = folder.split(".")[3:]

    # Now try opening the image, if it exists
    imagename = os.path.join(folder, os.path.basename(folder) +
                             ".clean.image.fits")
    feathered_imagename = \
        os.path.join(folder, os.path.basename(folder) +
                     ".clean.image.feathered.fits")

    if os.path.exists(feathered_imagename):
        # arr = fits.open(feathered_imagename)[0].data.squeeze()
        hi_hdu = fits.open(feathered_imagename)[0]
        interf_beam = Beam.from_fits_header(feathered_imagename)
    elif os.path.exists(imagename):
        # arr = fits.open(imagename)[0].data.squeeze()
        hi_hdu = fits.open(imagename)[0]
        interf_beam = Beam.from_fits_header(imagename)
    else:
        continue

    beam_factor = (sd_beam.sr / interf_beam.sr).value

    # Apply the clean mask.
    hi_hdu.data[~mask] = np.NaN

    feather_plot(hi_hdu, low_hdu, xaxisunit='arcsec')

    # fft_hi = np.fft.fftshift(np.fft.fft2(
    #     np.nan_to_num(arr * highresscalefactor)))

    # rad, azavg_hi_scaled = azimuthalAverage(np.abs(fft_hi * ikfft),
    #                                         returnradii=True)
    # rad, azavg_lo_scaled = azimuthalAverage(np.abs(fft_lo * kfft),
    #                                         returnradii=True)
    # rad, azavg_lo_deconv = azimuthalAverage(np.abs(fft_lo / kfft),
    #                                         returnradii=True)

    # # rad, azavg_hi_scaled = azimuthalAverage(np.abs(fft_hi),
    # #                                         returnradii=True)
    # # rad, azavg_lo_scaled = azimuthalAverage(np.abs(fft_lo),
    # #                                         returnradii=True)
    # # rad, azavg_lo_deconv = azimuthalAverage(np.abs(fft_lo),
    # #                                         returnradii=True)

    # OK = np.isfinite(azavg_kernel)

    # rad_pix = nax1 / rad
    # rad_as = pixscale * 3600 * rad_pix
    # if xaxisunit == 'lambda':
    #     lam = 1. / rad_as.to(u.rad).value
    #     xaxis = lam
    # elif xaxisunit == 'arcsec':
    #     xaxis = rad_as
    # else:
    #     raise ValueError("xaxisunit must be in (arcsec, lambda)")

    # ax1 = plt.subplot(1, 1, 1)
    # ax1.loglog(xaxis[OK], azavg_lo_scaled[OK], color='b', linewidth=2,
    #            alpha=0.5,
    #            linestyle='--',
    #            label="Low-res scaled image")
    # ax1.loglog(xaxis[OK], azavg_lo_deconv[OK], color='b', linewidth=2,
    #            alpha=0.5,
    #            linestyle=':',
    #            label="Low-res deconvolved image")
    # ax1.loglog(xaxis[OK], azavg_hi_scaled[OK], color='r', linewidth=2,
    #            alpha=0.5,
    #            linestyle='--',
    #            label="High-res scaled image")
    # plt.legend()

    plt.draw()
    print(split_name)
    raw_input("?")
    plt.clf()
