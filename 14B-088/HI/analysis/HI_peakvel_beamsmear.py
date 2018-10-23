
'''
Estimate the impact of beam smearing by running a standard deviation filter
over the peak velocity surface (since we use it for stacking) at different
resolutions.
'''

import numpy as np
import scipy.ndimage as nd
from astropy.io import fits
import astropy.units as u
from spectral_cube import Projection
import matplotlib.pyplot as plt
from turbustat.statistics import PDF
from scipy import stats
from scipy.stats import binned_statistic
import seaborn as sb
from corner import hist2d

from galaxy_params import gal_feath as gal

from paths import (fourteenB_wGBT_HI_file_dict, fourteenB_HI_data_wGBT_path,
                   allfigs_path)
from plotting_styles import onecolumn_figure, default_figure


def window_stdev(X, window_size):
    '''
    Standard deviation window.

    From: https://nickc1.github.io/python,/matlab/2016/05/17/Standard-Deviation-(Filters)-in-Matlab-and-Python.html
    '''
    r, c = X.shape
    X += np.random.rand(r, c) * 1e-6
    c1 = nd.uniform_filter(X, window_size, mode='reflect')
    c2 = nd.uniform_filter(X * X, window_size, mode='reflect')
    return np.sqrt(c2 - c1 * c1)


default_figure()

peakvels = \
    Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['PeakVels']))
mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0']))


# Set the window size to be the FWHM of the beam
# Pixels are essentially square in WCS here.

window_size = mom0.beam.major.to(u.deg) / \
    (peakvels.wcs.wcs.cdelt[-1] * u.deg)
window_size = int(np.ceil(window_size))

peakvels_val = peakvels.value.copy()
peakvels_val[np.isnan(peakvels_val)] = 0.0

stddev = window_stdev(peakvels_val, window_size)

# Now we need to mask out values near edges of NaNs
mask = np.isfinite(peakvels)
mask = nd.binary_erosion(mask, iterations=window_size)

# plt.imshow(stddev, origin='lower', vmax=1.e4)
# plt.contour(mask, colors='r')

onecolumn_figure()
sb.set_palette('colorblind')
col_pal = sb.color_palette('colorblind')

_ = plt.hist(np.log10(stddev[mask] / 1000.), bins='auto', alpha=0.4, label='19"')


peakvels_38 = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("smooth_2beam/M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.peakvels.fits")))
mom0_38 = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("smooth_2beam/M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.mom0.fits")))

# Set the window size to be the FWHM of the beam
# Pixels are essentially square in WCS here.

window_size_38 = mom0_38.beam.major.to(u.deg) / \
    (peakvels_38.wcs.wcs.cdelt[-1] * u.deg)
window_size_38 = int(np.ceil(window_size_38))

peakvels_val_38 = peakvels_38.value.copy()
peakvels_val_38[np.isnan(peakvels_val_38)] = 0.0

stddev_38 = window_stdev(peakvels_val_38, window_size_38)

# Now we need to mask out values near edges of NaNs
mask_38 = np.isfinite(peakvels_38)
mask_38 = nd.binary_erosion(mask_38, iterations=window_size_38)

mask_38[np.isnan(stddev_38)] = False
mask_38[stddev_38 == 0.0] = False

# Only compare where the original resolution values are defined
mask_38[~mask] = False

_ = plt.hist(np.log10(stddev_38[mask_38] / 1000.), bins='auto', alpha=0.4, label='38"')


peakvels_95 = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("smooth_5beam/M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.peakvels.fits")))
mom0_95 = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("smooth_5beam/M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.mom0.fits")))

# Set the window size to be the FWHM of the beam
# Pixels are essentially square in WCS here.

window_size_95 = mom0_95.beam.major.to(u.deg) / \
    (peakvels_95.wcs.wcs.cdelt[-1] * u.deg)
window_size_95 = int(np.ceil(window_size_95))

peakvels_val_95 = peakvels_95.value.copy()
peakvels_val_95[np.isnan(peakvels_val_95)] = 0.0

stddev_95 = window_stdev(peakvels_val_95, window_size_95)

# Now we need to mask out values near edges of NaNs
mask_95 = np.isfinite(peakvels_95)
mask_95 = nd.binary_erosion(mask_95, iterations=window_size_95)

mask_95[np.isnan(stddev_95)] = False
mask_95[stddev_95 == 0.0] = False

mask_95[~mask] = False

_ = plt.hist(np.log10(stddev_95[mask_95] / 1000.), bins='auto', alpha=0.4, label='95"')

plt.legend(frameon=True)
plt.xlim([-1., 2.])
plt.axvline(np.log10(0.2), linestyle='--', linewidth=3, alpha=0.6,
            color=col_pal[3])

# Fit some log-normals (close-ish to the shape) to get the peak location
# and estimate dispersion.

print("Fitting log-normals")

stddev_masked = stddev.copy()
stddev_masked[~mask] = np.NaN

# pdf = PDF(stddev_masked).run(verbose=False)

stddev_masked_38 = stddev_38.copy()
stddev_masked_38[~mask_38] = np.NaN

# pdf_38 = PDF(stddev_masked_38).run(verbose=False)

stddev_masked_95 = stddev_95.copy()
stddev_masked_95[~mask_95] = np.NaN

# pdf_95 = PDF(stddev_masked_95).run(verbose=False)


# def lognorm_mean(params):
#     return float(stats.lognorm.stats(s=params[0], scale=params[1],
#                  moments='m'))

# print("Mean: 19'' {0}; 38'' {1}; 95'' {2} km/s"
#       .format(lognorm_mean(pdf.model_params),
#               lognorm_mean(pdf_38.model_params),
#               lognorm_mean(pdf_95.model_params)))
# print("Dispersion: 19'' {0}; 38'' {1}; 95'' {2} km/s"
#       .format(pdf.model_params[1], pdf_38.model_params[1],
#               pdf_95.model_params[1]))

# plt.axvline(np.log10(lognorm_mean(pdf.model_params) / 1.e3), color=col_pal[0])
# plt.axvline(np.log10(lognorm_mean(pdf_38.model_params) / 1.e3), color=col_pal[1])
# plt.axvline(np.log10(lognorm_mean(pdf_95.model_params) / 1.e3), color=col_pal[2])

plt.axvline(np.log10(np.nanmedian(stddev_masked) / 1.e3), color=col_pal[0])
plt.axvline(np.log10(np.nanmedian(stddev_masked_38) / 1.e3), color=col_pal[1])
plt.axvline(np.log10(np.nanmedian(stddev_masked_95) / 1.e3), color=col_pal[2])

print("Median: 19'' {0}; 38'' {1}; 95'' {2} m/s"
      .format(np.nanmedian(stddev_masked),
              np.nanmedian(stddev_masked_38),
              np.nanmedian(stddev_masked_95)))
# Median: 19'' 1200.64344986; 38'' 1040.18519281; 95'' 2125.5324707 m/s

plt.grid()
plt.xlabel("log Standard Deviation (km/s)")
plt.tight_layout()
plt.savefig(allfigs_path('HI_properties/peakvel_stddevfilter_histograms.pdf'))
plt.savefig(allfigs_path('HI_properties/peakvel_stddevfilter_histograms.png'))
plt.close()

# Trends against each other?
# comb_mask_38 = np.logical_and(np.isfinite(stddev_masked),
#                               np.isfinite(stddev_masked_38))
# hist2d(stddev_masked[comb_mask_38], stddev_masked_38[comb_mask_38])

# comb_mask_95 = np.logical_and(np.isfinite(stddev_masked),
#                               np.isfinite(stddev_masked_95))
# hist2d(stddev_masked[comb_mask_95], stddev_masked_95[comb_mask_95])

# Against radius?
radius = gal.radius(header=mom0.header).to(u.kpc)

# hist2d(radius.value[comb_mask_38],
#        np.log10((stddev_masked_38 / stddev_masked)[comb_mask_38]))
# hist2d(radius.value[comb_mask_95],
#        np.log10((stddev_masked_95 / stddev_masked)[comb_mask_95]))

# Mild increase with radius for the inner 2 kpc.

# Create radial profiles of avg std

beam_pix = 41.

rad_bins = np.arange(0, 7.5, 0.5)
med_bin, bin_edges, cts = binned_statistic(radius.value[mask],
                                           stddev_masked[mask] / 1.e3,
                                           bins=rad_bins,
                                           statistic=np.mean)
# Last bin value are the points that are outside the largest bin
bin_cts = np.array([np.sum(cts == bin_lab) for bin_lab in np.unique(cts)[:-1]])

std_bin = binned_statistic(radius.value[mask],
                           stddev_masked[mask] / 1.e3,
                           bins=rad_bins,
                           statistic=np.std)[0]
# Correct by number of independent samples
# std_bin /= np.sqrt(bin_cts / beam_pix)

med_bin_38, bin_edges, cts = binned_statistic(radius.value[mask_38],
                                              stddev_masked_38[mask_38] / 1.e3,
                                              bins=rad_bins,
                                              statistic=np.mean)
bin_cts_38 = np.array([np.sum(cts == bin_lab)
                       for bin_lab in np.unique(cts)[:-1]])

std_bin_38 = binned_statistic(radius.value[mask_38],
                              stddev_masked_38[mask_38] / 1.e3,
                              bins=rad_bins,
                              statistic=np.std)[0]
# std_bin_38 /= np.sqrt(bin_cts_38 / beam_pix)

med_bin_95, bin_edges, cts = binned_statistic(radius.value[mask_95],
                                              stddev_masked_95[mask_95] / 1.e3,
                                              bins=rad_bins,
                                              statistic=np.mean)
bin_cts_95 = np.array([np.sum(cts == bin_lab)
                       for bin_lab in np.unique(cts)[:-1]])

std_bin_95 = binned_statistic(radius.value[mask_95],
                              stddev_masked_95[mask_95] / 1.e3,
                              bins=rad_bins,
                              statistic=np.std)[0]
# std_bin_95 /= np.sqrt(bin_cts_95 / beam_pix)

bin_cents = (bin_edges[1:] + bin_edges[:-1]) / 2.

plt.axhline(0.2, color=col_pal[3], linestyle='--', linewidth=4, alpha=0.75)
plt.axhline(2.6, color=col_pal[5], linestyle=':', linewidth=4, alpha=0.75)
plt.errorbar(bin_cents, med_bin, fmt='o-', drawstyle='steps-mid',
             yerr=std_bin, label='80 pc (19")')
plt.errorbar(bin_cents, med_bin_38, fmt='D-', drawstyle='steps-mid',
             yerr=std_bin_38, label='160 pc (38")')
plt.errorbar(bin_cents, med_bin_95, fmt='^-', drawstyle='steps-mid',
             yerr=std_bin_95, label='380 pc (95")')
plt.legend(frameon=True)
plt.grid()
plt.ylabel("Standard deviation of\n peak velocity (km/s)")
plt.xlabel("Radius (kpc)")
plt.tight_layout()
plt.savefig(allfigs_path('HI_properties/peakvel_stddevfilter_radprofile.pdf'))
plt.savefig(allfigs_path('HI_properties/peakvel_stddevfilter_radprofile.png'))
plt.close()


# Calculate the bin-averaged velocity dispersions from the radial profiles

pix_area = np.sum(radius < 7 * u.kpc)

med_bin_area_avg = (bin_cts * med_bin).sum() / pix_area
eight5_bin_area_avg = (bin_cts * (med_bin + std_bin)).sum() / pix_area
fifteen_bin_area_avg = (bin_cts * (med_bin - std_bin)).sum() / pix_area
med_bin_38_area_avg = (bin_cts * med_bin_38).sum() / pix_area
eight5_bin_38_area_avg = (bin_cts * (med_bin + std_bin_38)).sum() / pix_area
fifteen_bin_38_area_avg = (bin_cts * (med_bin - std_bin_38)).sum() / pix_area
med_bin_95_area_avg = (bin_cts * med_bin_95).sum() / pix_area
eight5_bin_95_area_avg = (bin_cts * (med_bin + std_bin_95)).sum() / pix_area
fifteen_bin_95_area_avg = (bin_cts * (med_bin - std_bin_95)).sum() / pix_area

print("Avg. broadening: {0}. Limits: {1}, {2}".format(med_bin_area_avg,
                                                      fifteen_bin_area_avg,
                                                      eight5_bin_area_avg))
print("38''Avg. broadening: {0}. Limits: {1}, {2}".format(med_bin_38_area_avg,
                                                          fifteen_bin_38_area_avg,
                                                          eight5_bin_38_area_avg))
print("95''Avg. broadening: {0}. Limits: {1}, {2}".format(med_bin_95_area_avg,
                                                          fifteen_bin_95_area_avg,
                                                          eight5_bin_95_area_avg))

# Can the line width increases at 95'' be entirely due to beam smearing?

hi_width_95 = 8.9
co_width_95 = 7.3

hi_width_38 = 8.0
co_width_38 = 6.0

hi_width_95_bc = np.sqrt(hi_width_95**2 - med_bin_95_area_avg**2)
hi_width_95_bc_err = (2 / hi_width_95_bc) * \
    (hi_width_95 * 0.1 + med_bin_95_area_avg * 1.0)

co_width_95_bc = np.sqrt(co_width_95**2 - med_bin_95_area_avg**2)
co_width_95_bc_err = (2 / co_width_95_bc) * \
    (co_width_95 * 1.3 + med_bin_95_area_avg * 1.0)

print("Beam smear-corr 95'' HI: {0}+/-{1}".format(hi_width_95_bc,
                                                  hi_width_95_bc_err))
print("Beam smear-corr 95'' CO: {0}+/-{1}".format(co_width_95_bc,
                                                  co_width_95_bc_err))
