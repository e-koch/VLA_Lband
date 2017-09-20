
'''
Figures from the galaxy-scale PV-slices made in HI_pvslices.py
'''

from spectral_cube import SpectralCube, Projection
import pvextractor as pv
from astropy.io import fits
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import mad_std
from aplpy import FITSFigure
import os

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path)
from galaxy_params import gal_feath as gal
from constants import distance, hi_freq
from plotting_styles import twocolumn_figure, default_figure

from HI_pvslices import phys_to_ang, obs_radius


cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.1kms.fits"))
mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0])
peakvels = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict["PeakVels"])[0])
mom1 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict["Moment1"])[0])
rotmod = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_noradial_nowarp_output/rad.fitmod.fits"))[0])

jybeam_to_K = cube.beam.jtok(hi_freq)

pa_angles = gal.position_angles(header=mom0.header)
radius = gal.radius(header=mom0.header)

# Set the angles to create PV slices from. Angles defined wrt to the disk PA.
# The cuts at  PA - 165 and 175 are to match the slices presented in Kam+17
thetas = np.array([0, 45, 90, 135, gal.position_angle.value - 165,
                   gal.position_angle.value - 175]) * u.deg
max_rad = 8.5 * u.kpc

pvslices = {}
paths = {}
for theta in thetas:

    # Adjust path length based on inclination
    obs_rad = obs_radius(max_rad, theta, gal)

    # Now convert the physical size in the observed frame to an angular size
    ang_length = 2 * phys_to_ang(obs_rad, distance)

    ang_width = 2 * phys_to_ang(obs_radius(max_rad, theta + 90 * u.deg, gal),
                                distance)

    pv_path = pv.PathFromCenter(gal.center_position, length=ang_length,
                                angle=theta + gal.position_angle,
                                sample=50,
                                width=ang_width)
    paths[int(theta.value)] = pv_path

    filename = "downsamp_1kms/M33_14B-088_HI.GBT_feathered_PA_{}_pvslice.fits".format(int(theta.value))

    try:
        proj = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path(filename))[0])
    except OSError:
        print("No {} saved.".format(filename))

    pvslices[int(theta.value)] = proj

rotsub_name = fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.GBT_feathered.rotation_corrected_PA_0_pvslice.fits")
pvslice_rotsub = Projection.from_hdu(fits.open(rotsub_name)[0])

# Make a plot of the slice along the major axis.

twocolumn_figure()

fig = plt.figure(figsize=(8.4, 4.02))

# major = pvslices[0]
major = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.GBT_feathered_PA_0_pvslice_200arcsec_width.fits"))[0])
major_path = paths[0]

major_std = mad_std(major)

fig1 = FITSFigure(major.hdu, figure=fig, subplot=(2, 1, 1))
fig1.show_grayscale(invert=True, stretch='arcsinh')  # , vmax=0.005)

contour_levels = major_std * [2, 4, 5]
fig1.show_contour(major.hdu,
                  levels=[1 / jybeam_to_K.value,
                          2 / jybeam_to_K.value,
                          3 / jybeam_to_K.value],
                  smooth=3)

# Levels in K
print("Contour levels in K: {}".format(contour_levels * jybeam_to_K))

# Now get the peak velocities from the path
xy_posns = \
    np.floor(np.round(np.array(major_path.get_xy(cube.wcs)))).astype('int')
y = xy_posns[:, 1]
x = xy_posns[:, 0]

peak_vel_pts = peakvels[y, x]
cent_vel_pts = mom1[y, x]
rot_vel_pts = rotmod[y, x]

# Now transform into distance along the path.
# Pixels are square here
pix_scale = np.abs(cube.header["CDELT2"]) * u.deg
dist = np.sqrt((y - y[0])**2 + (x - x[0])**2) * pix_scale

# fig1.show_markers(dist.value, peak_vel_pts.value, edgecolor='r',
#                   facecolor='r',
#                   marker='o', s=40)

# fig1.show_markers(dist.value, cent_vel_pts.value, edgecolor='b',
#                   facecolor='b',
#                   marker='D', s=40)

fig1.show_markers(dist.value, rot_vel_pts.value, edgecolor='g', facecolor='g',
                  marker='^', s=40)

fig2 = FITSFigure(pvslice_rotsub.hdu, figure=fig, subplot=(2, 1, 2))
fig2.show_grayscale(invert=True, stretch='arcsinh', vmax=0.005)

majorrotsub_std = mad_std(pvslice_rotsub)
contour_levels_rotsub = majorrotsub_std * [2, 4, 5]

print("Contour levels in K: {}".format(contour_levels_rotsub * jybeam_to_K))

fig2.show_contour(pvslice_rotsub.hdu,
                  levels=[1 / jybeam_to_K.value,
                          2 / jybeam_to_K.value,
                          3 / jybeam_to_K.value],
                  smooth=3)

# Add a line at a velocity of 0 (or Vsys for fig1) and at M33's center
zero_vel_posn = pvslice_rotsub.header['CRVAL2'] / \
    np.abs(pvslice_rotsub.header['CDELT2'])

fig1._ax1.axhline(zero_vel_posn, color='k', linestyle='-.',
                  linewidth=1, alpha=0.75)
fig2._ax1.axhline(zero_vel_posn, color='k', linestyle='-.',
                  linewidth=1, alpha=0.75)

# Add line at M33's center
# Must be in the center, since the pv path is defined wrt to the center.
fig2._ax1.axvline(pvslice_rotsub.shape[1] / 2, color='k', linestyle='-.',
                  linewidth=1, alpha=0.75)
fig1._ax1.axvline(major.shape[1] / 2, color='k', linestyle='-.',
                  linewidth=1, alpha=0.75)

fig_path = allfigs_path("pvslices")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

fig.savefig(os.path.join(fig_path, "M33_14B-088_major_axis_pvslices_wrotsub.png"))
fig.savefig(os.path.join(fig_path, "M33_14B-088_major_axis_pvslices_wrotsub.pdf"))
plt.close()


# Make another of just the unsubtracted data
fig = plt.figure(figsize=(8.4, 2.4))

major_K = fits.PrimaryHDU(major.value * jybeam_to_K.value, major.header)

fig1 = FITSFigure(major_K, figure=fig, subplot=(1, 1, 1))
fig1.show_grayscale(invert=True, stretch='arcsinh')  # , vmax=0.005)

contour_levels = major_std * [2, 4, 5]
fig1.show_contour(major_K,
                  levels=[1, 2, 3],
                  smooth=3)
fig1.show_markers(dist.value, rot_vel_pts.value, edgecolor='g', facecolor='g',
                  marker='^', s=40)
fig1._ax1.axhline(zero_vel_posn, color='k', linestyle='-.',
                  linewidth=1, alpha=0.75)
fig1._ax1.axvline(major.shape[1] / 2, color='k', linestyle='-.',
                  linewidth=1, alpha=0.75)

# Convert velocity to km /s
fig1._ax1.set_yticklabels(np.array([-300000, -250000, -200000,
                                    -150000, -100000]) / 1000)
fig1.set_axis_labels(ylabel='Velocity (km/s)', xlabel="Offset (deg)")

fig1.add_colorbar()
fig1.colorbar.set_location('top')
fig1.colorbar.set_label_properties(size=11)

plt.tight_layout()

fig.savefig(os.path.join(fig_path, "M33_14B-088_major_axis_pvslice.png"))
fig.savefig(os.path.join(fig_path, "M33_14B-088_major_axis_pvslice.pdf"))
plt.close()

default_figure()
