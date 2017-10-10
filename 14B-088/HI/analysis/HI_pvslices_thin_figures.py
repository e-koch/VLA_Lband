
'''
Make a figure of the thin pv-slices stacked on top of each other.
'''

from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy import units as u
import numpy as np
from glob import glob
import os
from os.path import join as osjoin
import matplotlib.pyplot as plt
from aplpy import FITSFigure

from paths import (fourteenB_HI_data_wGBT_path, allfigs_path,
                  fourteenB_wGBT_HI_file_dict)
from constants import hi_freq
from plotting_styles import twocolumn_figure, default_figure


# Make sure the figure directory exists
fig_path = allfigs_path("pvslices")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

pvslice_dir = fourteenB_HI_data_wGBT_path("downsamp_1kms/")

# I need the beam in the cube to convert to K
cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.1kms.fits"))

jybeam_to_K = cube.beam.jtok(hi_freq)

del cube

# Get all pv-slices
filenames = glob(osjoin(pvslice_dir, "M33_14B-088_HI.clean.image.GBT_feathered.1kms_PA_*_pvslice_40.0arcsec_width.fits"))

# The slices go from a PA of 0 to 175 in increments of 5 deg
ordered_filenames = []
pas = np.arange(0, 180, 5)

for pa in pas:
    for fname in filenames:
        if "PA_{}_".format(pa) in fname:
            ordered_filenames.append(fname)
            break

# Want to put on a common scale. Grab the max from all slices.
max_val = 0
for fname in ordered_filenames:
    hdu = fits.open(fname)[0]

    max_slice_val = np.nanmax(hdu.data)
    if max_slice_val > max_val:
        max_val = max_slice_val

# Split into figures of 6
for i in range(6):

    fig = plt.figure(figsize=(8.1, 11.))

    for j, n in enumerate(np.arange(6 * i, 6 * (i + 1))):

        hdu = fits.open(ordered_filenames[n])[0]

        fig_n = FITSFigure(hdu, subplot=(6, 1, j + 1), figure=fig)

        fig_n.show_grayscale(invert=True, stretch='arcsinh')

        fig_n.show_contour(hdu,
                           levels=[2 / jybeam_to_K.value,
                                   3 / jybeam_to_K.value,
                                   4 / jybeam_to_K.value],
                           smooth=3)

        zero_vel_posn = hdu.header['CRVAL2'] / \
            np.abs(hdu.header['CDELT2'])

        fig_n._ax1.axhline(zero_vel_posn * 1000., color='k', linestyle='-.',
                           linewidth=1, alpha=0.75)

        # Add line at M33's center
        # Must be in the center, since the pv path is defined wrt to the center
        fig_n._ax1.axvline(hdu.shape[1] / 2, color='k', linestyle='-.',
                           linewidth=1, alpha=0.75)

        fig_n._ax1.set_yticklabels(np.array([-300000, -250000, -200000,
                                            -150000, -100000]) / 1000)

        # fig_n.set_axis_labels(ylabel='Velocity (km/s)')
        fig_n.hide_axis_labels()
        fig_n.hide_ytick_labels()

        # Put the PA in the upper corner
        if i < 4:
            fig_n.add_label(0.81, 0.8, "{} deg".format(int(pas[n])),
                            relative=True, size=14,
                            bbox={"boxstyle": "square", "facecolor": "w"})
        else:
            fig_n.add_label(0.2, 0.8, "{} deg".format(int(pas[n])),
                            relative=True, size=14,
                            bbox={"boxstyle": "square", "facecolor": "w"})

        if j != 5:
            fig_n.hide_xaxis_label()
            fig_n.hide_xtick_labels()

        # if j == 0:
        #     fig_n.add_colorbar()
        #     fig_n.colorbar.set_location('top')
        #     fig_n.colorbar.set_label_properties(size=11)


    fig.savefig(osjoin(fig_path, "M33_14B-088_pvslices_40arcsec_{}.png".format(i)))
    fig.savefig(osjoin(fig_path, "M33_14B-088_pvslices_40arcsec_{}.pdf".format(i)))
    plt.close()


Now make a figure of all of the pv-slice paths on the zeroth moment map

mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0])
mom0.quicklook()
mom0.FITSFigure.show_regions(osjoin(pvslice_dir, "M33_14B-088_HI.clean.image.GBT_feathered.1kms_pvslice_40.0arcsec_width.reg"))
mom0.FITSFigure.save(osjoin(fig_path, "M33_14B-088_pvslices_40arcsec_paths.png"))
mom0.FITSFigure.save(osjoin(fig_path, "M33_14B-088_pvslices_40arcsec_paths.pdf"))


# Make a smaller figure for the paper. Include pv-slices along the major,
# minor and warped major (135 deg) axes


twocolumn_figure()

fig = plt.figure(figsize=(8.4, 4.2))

max_size = fits.open(ordered_filenames[0])[0].shape[1]
maj_header = fits.open(ordered_filenames[0])[0].header

for i, pa in zip(range(3), [0, 90, 135]):

    idx = np.where(pas == pa)[0]

    hdu = fits.open(ordered_filenames[idx])[0]

    # Convert to K
    hdu = fits.PrimaryHDU(hdu.data * jybeam_to_K.value, hdu.header)

    # Reverse the direction of the 135 slice
    if i == 2:
        hdu = fits.PrimaryHDU(hdu.data[::-1], hdu.header)

    if pa != 0:
        # Match the major axis slice length to make each the same shape
        padded_slice = np.zeros((hdu.shape[0], max_size)) * np.NaN
        # Velocity axes will match
        pad_size = (max_size - hdu.shape[1]) / 2
        if hdu.shape[1] % 2 == 0:
            padded_slice[:, pad_size:max_size - pad_size] = hdu.data
        else:
            padded_slice[:, pad_size:max_size - pad_size - 1] = hdu.data
        hdu = fits.PrimaryHDU(padded_slice, maj_header)

    fig_n = FITSFigure(hdu, subplot=(3, 1, i + 1), figure=fig)

    fig_n.show_grayscale(invert=True, stretch='arcsinh', vmin=0)

    fig_n.show_contour(hdu,
                       levels=[2,
                               3,
                               4],
                       smooth=3)

    # zero_vel_posn = hdu.header['CRVAL2'] / \
    #     np.abs(hdu.header['CDELT2'])

    # fig_n._ax1.axhline(zero_vel_posn * 1000., color='k', linestyle='-.',
    #                    linewidth=1, alpha=0.75)

    # Add line at M33's center
    # Must be in the center, since the pv path is defined wrt to the center
    fig_n._ax1.axvline(hdu.shape[1] / 2, color='k', linestyle='-.',
                       linewidth=1, alpha=0.75)

    fig_n._ax1.set_yticklabels(np.array([-300000, -250000, -200000,
                                        -150000, -100000]) / 1000)

    # fig_n.set_axis_labels(ylabel='Velocity (km/s)')
    # fig_n.hide_axis_labels()
    # fig_n.hide_ytick_labels()
    if i == 1:
        fig_n.set_axis_labels(ylabel='Velocity (km/s)')
    else:
        fig_n.axis_labels.hide_y()

    if i < 2:
        fig_n.axis_labels.hide_x()
        fig_n.tick_labels.hide_x()
    else:
        fig_n.set_axis_labels(xlabel='Offset (deg)')

    # Put the PA in the upper corner
    fig_n.add_label(0.9, 0.75, "{} deg".format(int(pa)),
                    relative=True, size=12,
                    bbox={"boxstyle": "square", "facecolor": "w"})

    fig_n.add_colorbar()
    fig_n.colorbar.set_ticks([0, 5, 10, 20, 40])
    fig_n.colorbar.set_font(size=11)

    if i == 1:
        fig_n.colorbar.set_axis_label_text('Intensity (K)')
        fig_n.colorbar.set_axis_label_font(size=12)

plt.subplots_adjust(hspace=0.01)
plt.tight_layout()

fig.savefig(osjoin(fig_path, "M33_14B-088_pvslices_40arcsec_PA_0_90_135.png"))
fig.savefig(osjoin(fig_path, "M33_14B-088_pvslices_40arcsec_PA_0_90_135.pdf"))
plt.close()

default_figure()
