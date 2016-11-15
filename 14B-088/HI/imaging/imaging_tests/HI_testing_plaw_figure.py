
'''
Create a figure summarizing differences in the power spectra.

Originally setup to have 4 panels, but the multiscale and mask comparisons
are nearly identical. So only the w/ and w/o model and kernel weights are
returned.
'''

import os
from glob import glob
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from radio_beam import Beam
from reproject import reproject_interp
from FITS_tools.hcongrid import hcongrid
from scipy import ndimage as nd

from astropy.utils.console import ProgressBar
import astropy.io.fits as fits
from astropy import wcs
from image_tools.radialprofile import azimuthalAverage

data_path = os.path.expanduser("~/MyRAID/M33/14B-088/HI/channel_testing")
# data_path = os.path.expanduser("/Volumes/Travel_Data/M33/testing")


def append_path(path):
    return os.path.join(data_path, path)


def avg_fwhm(beam):
    return np.sqrt(beam.major.to(u.arcsec) * beam.minor.to(u.arcsec))


parameters = ["CASAVer", "Model", "Mask", "AllFields", "MScale", "Tclean"]

# Load in the mask
mask_hdu = fits.open(append_path("M33_14B-088_HI_mask_channel_330.fits"))[0]
mask = mask_hdu.data.squeeze()
mask = mask[::-1, ::-1] == 1

low_hdu = fits.open(append_path("M33_14B-088_HI_model_channel_330.fits"))[0]
model = low_hdu.data.squeeze()
model = model[::-1, ::-1]

sd_beam = Beam.from_fits_header(append_path("M33_14B-088_HI_model_channel_330.fits"))
fft_sd = np.fft.fftshift(np.fft.fft2(np.nan_to_num(model)))

# Also load in the original Arecibo channel that has not been regridded.
low_hdu_orig = \
    fits.open(append_path("M33_14B-088_HI_model_original_arecibo.fits"))[0]

# The commented out sections were for ensuring that the model matched
# regridded versions using FITS_tools and reproject. On the large scales that
# matter, the answer is yes they do match.
# regrid_model, footprint = reproject_interp(low_hdu_orig, low_hdu.header)
# regrid_model_fitstools = hcongrid(low_hdu_orig.data, low_hdu_orig.header,
#                                   low_hdu.header)
regrid_mask = reproject_interp(mask_hdu, low_hdu_orig.header)[0]
regrid_mask[np.isnan(regrid_mask)] = 0.0
regrid_mask = regrid_mask == 1.0
regrid_mask = nd.binary_fill_holes(regrid_mask)


low_hdu_orig.data[~regrid_mask] = np.NaN
# regrid_model[~mask[::-1, ::-1]] = np.NaN

fft_sd_arecibo = np.fft.fftshift(np.fft.fft2(np.nan_to_num(low_hdu_orig.data)))

pixscale_arecibo = \
    wcs.utils.proj_plane_pixel_area(wcs.WCS(low_hdu_orig))**0.5 * \
    u.deg

# SD Azimuthal average
rad, azavg_sd = azimuthalAverage(np.abs(fft_sd),
                                 returnradii=True)

rad_arecibo, azavg_sd_arecibo = \
    azimuthalAverage(np.abs(fft_sd_arecibo), returnradii=True)
# rad_model, azavg_sd_model = \
#     azimuthalAverage(np.abs(np.fft.fftshift(np.fft.fft2(np.nan_to_num(regrid_model)))), returnradii=True)
# azavg_sd_model_fitstools = \
#     azimuthalAverage(np.abs(np.fft.fftshift(np.fft.fft2(np.nan_to_num(regrid_model_fitstools)))))

# Grid scaling used.
pixscale = 3 * u.arcsec
# Image shape
nax2, nax1 = model.shape

# High scale factor
highresscalefactor = 1.0
# Low scale factor
lowresscalefactor = 1.0

sd_kernel = sd_beam.as_kernel(pixscale, x_size=nax2, y_size=nax1)
kfft = np.abs(np.fft.fftshift(np.fft.fft2(sd_kernel)))

azavg_kernel = azimuthalAverage(np.abs(kfft))

OK = np.isfinite(azavg_kernel)

rad_pix = nax1 / rad
xaxis = pixscale.to(u.arcsec).value * rad_pix

rad_pix_arecibo = max(fft_sd_arecibo.shape) / rad_arecibo
xaxis_arecibo = pixscale_arecibo.to(u.arcsec).value * rad_pix_arecibo

# w/ and w/o model

prefix = "14B-088_HI_LSRK.ms.contsub_channel_1000.CASAVer_440.Model_{0}" \
    ".Mask_{1}.AllFields_{2}.MScale_{3}.Tclean_F"

wo_model_name = prefix.format('F', 'T', 'T', 'T')
wo_model_hdu = fits.open(append_path("{0}/{0}.clean.image.fits".format(wo_model_name)))[0]
w_model_name = prefix.format('T', 'T', 'T', 'T')
w_model_hdu = \
    fits.open(append_path("{0}/{0}.clean.image.feathered.fits".format(w_model_name)))[0]

interf_beam = Beam.from_fits_header(wo_model_hdu.header)
beam_factor = (sd_beam.sr / interf_beam.sr).value

# And the factor to the original Arecibo grid
beam_factor_arecibo = beam_factor * \
    (pixscale / pixscale_arecibo.to(u.arcsec)).value ** 2

# Calculate the FWHMs
avg_sd_fwhm = avg_fwhm(sd_beam).value
avg_interf_fwhm = avg_fwhm(interf_beam).value

# The shortest baseline in the data is 44.777 m
short_base = 44.777 * u.m
lambda_hi = (1420.40575177 * u.MHz).to(u.m, equivalencies=u.spectral())
short_base_scale = (lambda_hi / short_base)\
    .to(u.arcsec, equivalencies=u.dimensionless_angles()).value

azavg_sd /= beam_factor
azavg_sd_arecibo /= beam_factor_arecibo

# Apply the clean mask.
wo_model_hdu.data[~mask] = np.NaN
w_model_hdu.data[~mask] = np.NaN

fft_wo_model = np.fft.fftshift(np.fft.fft2(np.nan_to_num(wo_model_hdu.data)))
fft_w_model = np.fft.fftshift(np.fft.fft2(np.nan_to_num(w_model_hdu.data)))

azavg_wo_model = azimuthalAverage(np.abs(fft_wo_model))
azavg_w_model = azimuthalAverage(np.abs(fft_w_model))

# ax1 = plt.subplot(2, 2, 1)
ax1 = plt.subplot(2, 2, 1)
ax1.loglog(xaxis[OK], azavg_wo_model[OK], color='b', linewidth=3,
           alpha=0.5,
           linestyle='-.',
           label="VLA")
ax1.loglog(xaxis[OK], azavg_sd[OK], color='r', linewidth=2,
           alpha=0.5,
           linestyle='--',
           label="Arecibo")
ax1.loglog(xaxis_arecibo, azavg_sd_arecibo, color='g', linewidth=2,
           alpha=0.5,
           linestyle='--',
           label="Arecibo Orig")
ax1.loglog(xaxis[OK], azavg_w_model[OK], color='k', linewidth=2,
           alpha=1.0,
           linestyle='-',
           label="VLA + Arecibo")
# ax1.set_xlabel("Angular scale (arcsec)")
ax1.set_ylabel("Power spectrum $|FT|$")

ax1.grid(True)

# Cut off tiny scale beyond VLA resolution
ax1.set_xlim([10, 1.2e4])

# Set by-eye
ylims = [1e-4, 2.5e3]
ax1.set_ylim(ylims)

# Add vertical lines at the beam's major FWHMs
ax1.vlines(avg_sd_fwhm, ylims[0], ylims[1],
           colors='gray')
ax1.vlines(avg_interf_fwhm, ylims[0], ylims[1],
           colors='gray')
ax1.vlines(short_base_scale, ylims[0], ylims[1],
           colors='gray', linestyle='--')

ax1.legend(loc='lower right', frameon=True)


# With dirty VLA only

wo_model_hdu = fits.open(append_path("14B-088_HI_LSRK.ms.contsub_channel_1000.dirty_imaging/"
                          "14B-088_HI_LSRK.ms.contsub_channel_1000.clean.image.fits"))[0]
w_model_name = prefix.format('T', 'T', 'T', 'T')
w_model_hdu = \
    fits.open(append_path("{0}/{0}.clean.image.feathered.fits".format(w_model_name)))[0]

# Apply the clean mask.
wo_model_hdu.data[~mask] = np.NaN
w_model_hdu.data[~mask] = np.NaN

fft_wo_model = np.fft.fftshift(np.fft.fft2(np.nan_to_num(wo_model_hdu.data)))
fft_w_model = np.fft.fftshift(np.fft.fft2(np.nan_to_num(w_model_hdu.data)))

azavg_wo_model = azimuthalAverage(np.abs(fft_wo_model))
azavg_w_model = azimuthalAverage(np.abs(fft_w_model))

# ax1 = plt.subplot(2, 2, 1)
ax2 = plt.subplot(2, 2, 2)
ax2.loglog(xaxis[OK], azavg_wo_model[OK], color='b', linewidth=3,
           alpha=0.5,
           linestyle='-.',
           label="VLA Dirty")
ax2.loglog(xaxis[OK], azavg_sd[OK], color='r', linewidth=2,
           alpha=0.5,
           linestyle='--',
           label="Arecibo")
ax2.loglog(xaxis[OK], azavg_w_model[OK], color='k', linewidth=2,
           alpha=1.0,
           linestyle='-',
           label="VLA + Arecibo")
# ax2.set_xlabel("Angular scale (arcsec)")
# ax2.set_ylabel("Power spectrum $|FT|$")

ax2.grid(True)

# Cut off tiny scale beyond VLA resolution
ax2.set_xlim([10, 1.2e4])

# Set by-eye
ylims = [1e-4, 2.5e3]
ax2.set_ylim(ylims)

# Add vertical lines at the beam's major FWHMs
ax2.vlines(avg_sd_fwhm, ylims[0], ylims[1],
           colors='gray')
ax2.vlines(avg_interf_fwhm, ylims[0], ylims[1],
           colors='gray')
ax2.vlines(short_base_scale, ylims[0], ylims[1],
           colors='gray', linestyle='--')

ax2.legend(loc='lower right', frameon=True)

# Next plot is w/ model, w/ and w/o multiscale

# wo_mscale_name = prefix.format('T', 'T', 'T', 'F')
# wo_mscale_hdu = fits.open(append_path("{0}/{0}.clean.image.feathered.fits".format(wo_mscale_name)))[0]
# w_mscale_name = prefix.format('T', 'T', 'T', 'T')
# w_mscale_hdu = \
#     fits.open(append_path("{0}/{0}.clean.image.feathered.fits".format(w_mscale_name)))[0]

# # Apply the clean mask.
# wo_mscale_hdu.data[~mask] = np.NaN
# w_mscale_hdu.data[~mask] = np.NaN

# fft_wo_mscale = np.fft.fftshift(np.fft.fft2(np.nan_to_num(wo_mscale_hdu.data)))
# fft_w_mscale = np.fft.fftshift(np.fft.fft2(np.nan_to_num(w_mscale_hdu.data)))

# azavg_wo_mscale = azimuthalAverage(np.abs(fft_wo_mscale))
# azavg_w_mscale = azimuthalAverage(np.abs(fft_w_mscale))

# ax2 = plt.subplot(2, 2, 2)
# ax2.loglog(xaxis[OK], azavg_wo_mscale[OK], color='b', linewidth=3,
#            alpha=0.5,
#            linestyle='-.',
#            label="Without multiscale")
# ax2.loglog(xaxis[OK], azavg_w_mscale[OK], color='r', linewidth=2,
#            alpha=0.5,
#            linestyle='--',
#            label="With multiscale")

# ax2.grid(True)

# # Cut off tiny scale beyond VLA resolution
# ax2.set_xlim([10, 1.2e4])

# # Set by-eye
# ylims = [7e-2, 2.5e3]
# ax2.set_ylim(ylims)

# # Add vertical lines at the beam's major FWHMs
# ax2.vlines(avg_sd_fwhm, ylims[0], ylims[1],
#            colors='gray')
# ax2.vlines(avg_interf_fwhm, ylims[0], ylims[1],
#            colors='gray')

# ax2.legend(loc='lower right', frameon=True)

# Kernel weightings

inter_kernel = interf_beam.as_kernel(pixscale, x_size=nax2, y_size=nax1)
kfft_inter = np.abs(np.fft.fftshift(np.fft.fft2(inter_kernel)))
ikfft = 1 - kfft_inter

azavg_kernel = azimuthalAverage(np.abs(kfft))
azavg_ikernel = azimuthalAverage(np.abs(ikfft))

# ax3 = plt.subplot(2, 2, 3)
ax3 = plt.subplot(2, 1, 2)
ax3.semilogx(xaxis[OK], azavg_kernel[OK], color='b', linewidth=2, alpha=0.8,
           label="Arecibo beam")
ax3.semilogx(xaxis[OK], azavg_ikernel[OK], color='r', linewidth=2, alpha=0.8,
           label="VLA beam", linestyle='--')

ylims = [1e-5, 1.1]
ax3.vlines(avg_sd_fwhm, ylims[0], ylims[1],
           colors='gray')
ax3.vlines(avg_interf_fwhm, ylims[0], ylims[1],
           colors='gray')
ax3.vlines(short_base_scale, ylims[0], ylims[1],
           colors='gray', linestyle='--')

ax3.set_ylim([0,1])

ax3.set_xlim([10, 1.2e4])

ax3.set_ylabel("Kernel Weight")
ax3.set_xlabel("Angular scale (arcsec)")

ax3.grid(True)
ax3.legend(loc='lower right', frameon=True)

# Next plot is w/ model, w/ and w/o mask

# wo_mask_name = prefix.format('T', 'F', 'T', 'T')
# wo_mask_hdu = fits.open(append_path("{0}/{0}.clean.image.feathered.fits".format(wo_mask_name)))[0]
# w_mask_name = prefix.format('T', 'T', 'T', 'T')
# w_mask_hdu = \
#     fits.open(append_path("{0}/{0}.clean.image.feathered.fits".format(w_mask_name)))[0]

# # Apply the clean mask.
# wo_mask_hdu.data[~mask] = np.NaN
# w_mask_hdu.data[~mask] = np.NaN

# fft_wo_mask = np.fft.fftshift(np.fft.fft2(np.nan_to_num(wo_mask_hdu.data)))
# fft_w_mask = np.fft.fftshift(np.fft.fft2(np.nan_to_num(w_mask_hdu.data)))

# azavg_wo_mask = azimuthalAverage(np.abs(fft_wo_mask))
# azavg_w_mask = azimuthalAverage(np.abs(fft_w_mask))

# ax4 = plt.subplot(2, 2, 4)
# ax4.loglog(xaxis[OK], azavg_wo_mask[OK], color='b', linewidth=3,
#            alpha=0.5,
#            linestyle='-.',
#            label="Without mask")
# ax4.loglog(xaxis[OK], azavg_w_mask[OK], color='r', linewidth=2,
#            alpha=0.5,
#            linestyle='--',
#            label="With mask")
# ax4.set_xlabel("Angular scale (arcsec)")

# ax4.grid(True)

# # Cut off tiny scale beyond VLA resolution
# ax4.set_xlim([10, 1.2e4])

# # Set by-eye
# ylims = [7e-2, 2.5e3]
# ax4.set_ylim(ylims)

# # Add vertical lines at the beam's major FWHMs
# ax4.vlines(avg_sd_fwhm, ylims[0], ylims[1],
#            colors='gray')
# ax4.vlines(avg_interf_fwhm, ylims[0], ylims[1],
#            colors='gray')

# ax4.legend(loc='lower right', frameon=True)

# ax4.set_xlabel("Angular scale (arcsec)")
