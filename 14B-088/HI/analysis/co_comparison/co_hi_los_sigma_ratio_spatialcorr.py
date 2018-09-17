
'''
Look for typical correlation scales in the line width ratio map.
'''


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.utils.console import ProgressBar
import os
from itertools import product

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path,
                   iram_co21_14B088_data_path)
from plotting_styles import default_figure, onecolumn_twopanel_figure
from galaxy_params import gal_feath as gal

default_figure()

fig_path = allfigs_path("co_vs_hi")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_co_gaussfit_column_densities_perpix.fits"))

# Don't consider the "bad fits" that are probably due to multiple components
good_pts = np.logical_and(~tab['multicomp_flag_HI'],
                          ~tab['multicomp_flag_CO'])
good_pts = np.logical_and(good_pts,
                          tab["sigma_HI"] > 3800)
# Minimum CO line width of one channel.
good_pts = np.logical_and(good_pts,
                          tab["sigma_CO"] >= 2600)


moment0 = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]

ratio_map = np.zeros(moment0.data.shape)

ratio_map[tab['ypts'][good_pts], tab['xpts'][good_pts]] = \
    (tab['sigma_CO'] / tab['sigma_HI'])[good_pts]

ratio_map[ratio_map == 0.] = np.NaN

def twopt_corr(image, max_radius=51, boundary='cut'):

    pix_lags = (np.arange(max_radius) - max_radius // 2)

    dx = pix_lags.copy()
    dy = pix_lags.copy()

    corr_surface = np.zeros((max_radius, max_radius))

    bar = ProgressBar(len(dx) * len(dy))

    for n, (x_shift, y_shift) in enumerate(product(dx, dy)):

        i, j = np.unravel_index(n, (len(dx), len(dy)))

        if x_shift == 0 and y_shift == 0:
            corr_surface[j, i] = 1.

        if x_shift == 0:
            tmp = image
        else:
            tmp = np.roll(image, x_shift, axis=0)

        if y_shift != 0:
            tmp = np.roll(tmp, y_shift, axis=1)

        if boundary is "cut":
            # Always round up to the nearest integer.
            x_shift = np.ceil(x_shift).astype(int)
            y_shift = np.ceil(y_shift).astype(int)
            if x_shift < 0:
                x_slice_data = slice(None, tmp.shape[0] + x_shift)
                x_slice_tmp = slice(-x_shift, None)
            else:
                x_slice_data = slice(x_shift, None)
                x_slice_tmp = slice(None, tmp.shape[0] - x_shift)

            if y_shift < 0:
                y_slice_data = slice(None, tmp.shape[1] + y_shift)
                y_slice_tmp = slice(-y_shift, None)
            else:
                y_slice_data = slice(y_shift, None)
                y_slice_tmp = slice(None, tmp.shape[1] - y_shift)

            data_slice = (x_slice_data, y_slice_data)
            tmp_slice = (x_slice_tmp, y_slice_tmp)
        elif boundary is "continuous":
            data_slice = (slice(None),) * 2
            tmp_slice = (slice(None),) * 2

        # Remove NaNs
        img_vals = image[data_slice]
        tmp_vals = tmp[tmp_slice]

        good_posns = np.logical_and(np.isfinite(img_vals),
                                    np.isfinite(tmp_vals))

        corr_surface[j, i] = np.corrcoef(img_vals[good_posns],
                                         tmp_vals[good_posns])[0, 1]

        bar.update(n + 1)

    return corr_surface


ratio_corr = twopt_corr(ratio_map)

# Estimate the width from the 1 / e contour level
from turbustat.statistics.pca.width_estimate import \
    (get_contour_path, fit_2D_ellipse)

max_radius = 51
pix_lags = (np.arange(max_radius) - max_radius // 2)

ymat, xmat = np.meshgrid(pix_lags, pix_lags, indexing='ij')

level = np.exp(-1)
paths = get_contour_path(xmat, ymat, ratio_corr, level)
pidx = np.where([p.contains_point((0, 0)) for p in paths])[0]

good_path = paths[pidx[0]]

output = fit_2D_ellipse(good_path.vertices)

scales = np.sqrt(output[0]**2 + output[1]**2)

scale_errors = \
    np.sqrt((output[0] * output[2])**2 +
            (output[1] * output[3])**2) / scales

print("Correlation scale: {0:.2f}+/-{1:.2f}".format(scales, scale_errors))
# Correlation scale: 2.19+/-0.02

# The line width ratios are correlated primarily on beam scales

plt.imshow(ratio_corr, vmin=0., vmax=1., cmap='viridis')
plt.colorbar()
