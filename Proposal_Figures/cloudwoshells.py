
'''
Figures and estimates for the severely HI poor molecular cloud.
'''

import numpy as np
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube import BooleanArrayMask
import matplotlib.pyplot as p
import astropy.io.fits as fits
import astropy.units as u
from astropy.table import Table
import os
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.wcs import WCS
import astropy.constants as const
from signal_id import Noise
from radio_beam import Beam
from astropy.utils.console import ProgressBar
from corner import hist2d
from aplpy import FITSFigure
from reproject import reproject_interp

from paths import c_hi_analysispath

os.sys.path.insert(0, c_hi_analysispath(""))

from paths import (root, fourteenB_HI_data_path, iram_co21_data_path,
                            proposal_figures_path)
from constants import (hi_freq, cube_name, hi_mass_conversion,
                       co21_mass_conversion, ang_to_phys,
                       moment0_name)
from galaxy_params import gal
from plotting_styles import (onecolumn_figure, default_figure,
                             twocolumn_figure, align_yaxis)
from co_comparison.co_cube_masking import simple_masking
from co_comparison.krumholz_2009_model import krumholz_ratio_model
from spectra_shifter import spectrum_shifter


def spatial_reproject(cube, header):

    celestial_header = WCS(header).celestial.to_header()
    celestial_header["NAXIS1"] = header["NAXIS1"]
    celestial_header["NAXIS2"] = header["NAXIS2"]
    celestial_header["NAXIS"] = 2

    newdata = np.empty((cube.shape[0], header["NAXIS2"], header["NAXIS1"]))
    for i in ProgressBar(cube.shape[0]):
        newdata[i] = cube[i].reproject(celestial_header).value

    newheader = cube.header.copy()
    for key in ["CDELT", "CRPIX", "CUNIT", "CTYPE", "CRVAL", "NAXIS"]:
        try:
            newheader[key + "1"] = header[key + "1"]
            newheader[key + "2"] = header[key + "2"]
        except KeyError:
            continue

    return SpectralCube(newdata, WCS(newheader), meta=cube.meta,
                        header=newheader).with_mask(np.isfinite(newdata))


hi_cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))
cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))

cloud_mask = fits.open(iram_co21_data_path('asgn.fits'))[0]
# Cloud 139's footprint
footprint = (cloud_mask.data == 139).sum(0) > 0
num_pix = footprint.sum()

# Load in the txt files.
cloud_avg_specs_co = np.loadtxt(fourteenB_HI_data_path("GMC_stackedspectra_peakvels_CO.txt"))
cloud_avg_specs_hi = np.loadtxt(fourteenB_HI_data_path("GMC_stackedspectra_peakvels_HI.txt"))

# In order to easily calculate moments, let's re-arrange the spectra into a
# cube
# Use the WCS from the actual cubes, since all we care about here is the
# spectral dimension.
co_shape = cloud_avg_specs_co.shape
co_spectra = SpectralCube(cloud_avg_specs_co.T[:, :, np.newaxis], wcs=cube.wcs,
                          meta={"BUNIT": "K"})

# Get the angular pixel size
ang_size = proj_plane_pixel_scales(cube.wcs)[0] * u.deg
phys_size = ang_to_phys(ang_size)
# The physical size of the cloud is then (at least the part used in the stack):
cloud_area = num_pix * phys_size**2

hi_shape = cloud_avg_specs_hi.shape
hi_spectra = SpectralCube(cloud_avg_specs_hi.T[:, :, np.newaxis],
                          wcs=hi_cube.wcs, meta={"BUNIT": "K"})

# Cut out the relevant velocity range (based on the HI)
co_spectra_cloud = co_spectra.spectral_slab(-173200 * u.m/u.s, -192300 * u.m / u.s)
hi_spectra_cloud = hi_spectra.spectral_slab(-173200 * u.m/u.s, -192300 * u.m / u.s)

co_mom0_clouds = co_spectra_cloud.moment0().to(u.K * u.km / u.s)
hi_mom0_clouds = hi_spectra_cloud.moment0().to(u.K * u.km / u.s)

# We're interested in cloud 139
cloud_ind = 138
co_cloud_mass = co_mom0_clouds[cloud_ind, 0] * co21_mass_conversion * cloud_area
hi_cloud_mass = hi_mom0_clouds[cloud_ind, 0] * hi_mass_conversion * cloud_area

print("CO Cloud mass {}".format(co_cloud_mass))
print("HI Cloud mass {}".format(hi_cloud_mass))

# To follow the rest of the population, need Sigma_T = 16 Msol/pc^2.

# The average H2-to-HI ratio of the cloud is ~0.25
tab = Table.read(fourteenB_HI_data_path("tables/column_densities_perpix.fits"))

# Construct the x-y points
low_pts = np.logical_and(tab["Sigma_Total"] < 4.5,
                         np.log10(tab["Ratio"]) > -0.15)
low_pts = np.logical_or(low_pts, tab["Sigma_Total"] < 1.8)

# There's one pixel from much farther south
low_pts = np.logical_and(low_pts, tab['xpix'] > 760)

# Separate the point belonging to cloud 139
cloud_pts = np.logical_and(low_pts, tab['xpix'] < 780)

total_sd_pix = tab["Sigma_Total"].quantity
gas_ratio_pix = tab["Ratio"].quantity
atomic_sd_pix = tab["Sigma_HI"].quantity
molec_sd_pix = tab["Sigma_H2"].quantity

avg_ratio = np.mean(gas_ratio_pix[low_pts])
avg_sigma = np.mean(total_sd_pix[low_pts])
avg_sigma_HI = np.mean(atomic_sd_pix[low_pts])
avg_sigma_H2 = np.mean(molec_sd_pix[low_pts])

# How much HI are we missing? And how many ionizing photons are needed to remove it?
N_HI = ((10 * u.solMass / u.pc**2 - avg_sigma_HI) * cloud_area) / \
    (1.4 * const.m_p.to(u.solMass))

# Overplot points on various images:
co_mom0_hdu = fits.open(iram_co21_data_path("m33.co21_iram.masked.moment0.hireprojection.fits"))[0]
co_mom0 = Projection.from_hdu(co_mom0_hdu)

default_figure()

twocolumn_figure()
fig = p.figure(figsize=(8.4, 6.65))

reg_slice = [slice(750, 800), slice(600, 650)]

props = dict(boxstyle='round', facecolor='white', alpha=1.0)

fig1 = FITSFigure(co_mom0.hdu)#[reg_slice].hdu)#, subplot=(2, 2, 1), figure=fig)
fig1.show_grayscale()
fig1.set_nan_color("black")
fig1.show_markers(tab["RA"][low_pts], tab["Dec"][low_pts],
                  edgecolor='blue', facecolor='blue',
                  marker='D', s=20, alpha=0.7)
fig1.show_markers(tab["RA"][posns_141], tab["Dec"][posns_141],
                  edgecolor='red', facecolor='red',
                  marker='D', s=20, alpha=0.7)
fig1.add_label(0.1, 0.1, "CO(2-1)", relative=True, bbox=props,
               horizontalalignment='left')
fig1.show_regions("hipoor_GMCs.reg")
fig1.hide_axis_labels()
fig1.hide_tick_labels()

hi_mom0_hdu = fits.open(fourteenB_HI_data_path(moment0_name))[0]
hi_mom0 = Projection.from_hdu(hi_mom0_hdu)
fig2 = FITSFigure(hi_mom0[reg_slice].hdu, subplot=(2, 2, 2),
                  figure=fig)
fig2.show_grayscale()
fig2.set_nan_color("black")
fig2.add_label(0.1, 0.1, "HI", relative=True, bbox=props,
               horizontalalignment='left')
fig2.show_regions("hipoor_GMCs.reg")
fig2.hide_axis_labels()
fig2.hide_tick_labels()

# Halpha w/ FUV contours
ha_hdu = fits.open("/home/eric/MyRAID/M33/ha/ha.fits")[0]
ha = Projection.from_hdu(ha_hdu)
lower_corner = hi_mom0.wcs.all_pix2world(600, 750, 0)
lower_corner_ha = ha.wcs.all_world2pix(lower_corner[0], lower_corner[1], 0)
upper_corner = hi_mom0.wcs.all_pix2world(650, 800, 0)
upper_corner_ha = ha.wcs.all_world2pix(upper_corner[0], upper_corner[1], 0)
reg_slice_ha = [slice(int(lower_corner_ha[1]), int(upper_corner_ha[1])),
                slice(int(lower_corner_ha[0]), int(upper_corner_ha[0]))]

fig3 = FITSFigure(ha[reg_slice_ha].hdu, subplot=(2, 2, 3), figure=fig)
fig3.show_grayscale(vmin=0, vmax=800)
fig3.show_contour("/home/eric/MyRAID/M33/Spitzer/mips1/mips24.fits",
                  levels=[2, 4.5, 15])
fig3.show_regions("hipoor_GMCs_band3.reg")
fig3.show_regions("hipoor_GMCs_xraysrcs.reg")
fig3.add_label(0.1, 0.1, r"H$\alpha$ with 24 $\mu$m contour", relative=True,
               bbox=props, horizontalalignment='left')
fig3.hide_axis_labels()
fig3.hide_tick_labels()

# B band
bband_hdu = fits.open("/home/eric/MyRAID/M33/Bband/M33CBm.qzap.fits.gz")[0]
bband = Projection.from_hdu(bband_hdu)
lower_corner_bband = bband.wcs.all_world2pix(lower_corner[0], lower_corner[1], 0)
upper_corner_bband = bband.wcs.all_world2pix(upper_corner[0], upper_corner[1], 0)
reg_slice_bb = [slice(int(lower_corner_bband[1]), int(upper_corner_bband[1])),
                slice(int(lower_corner_bband[0]), int(upper_corner_bband[0]))]

fig4 = FITSFigure(bband[reg_slice_bb].hdu, subplot=(2, 2, 4),
                  figure=fig)
fig4.show_grayscale()
fig4.show_contour("/home/eric/MyRAID/M33/fuv/galex_fuv.fits",
                  levels=[0.05, 0.1, 0.4, 0.7])
fig4.add_label(0.1, 0.1, r"B-band with FUV contour", relative=True,
               bbox=props, horizontalalignment='left')
fig4.hide_axis_labels()
fig4.hide_tick_labels()
fig4.show_regions("hipoor_GMCs_band3.reg")
fig4.show_regions("hipoor_GMCs_xraysrcs.reg")

p.tight_layout()
save_name = "hipoorGMCs_multiwave_almacoverage"
p.savefig(proposal_figures_path("{}.pdf".format(save_name)))
p.savefig(proposal_figures_path("{}.png".format(save_name)))
p.close()

# Remake the contour plot, and put the outlier points in another colour.
onecolumn_figure()

hist2d(total_sd_pix.value, np.log10(gas_ratio_pix.value),
       data_kwargs={"alpha": 0.6})
p.xlim([0, 60])
p.ylim([-2.1, 0.7])
p.xlabel("$\Sigma_{\mathrm{Gas}}$ (M$_{\odot}$ pc$^{-2}$)")
p.ylabel("log H$_2$-to-HI Ratio $\Sigma_{\mathrm{H2}} /"
         " \Sigma_{\mathrm{HI}}$")
p.plot(total_sd_pix[low_pts].value, np.log10(gas_ratio_pix[low_pts].value),
       'bD', ms=3)

sds = np.arange(0.1, 60, 0.2)
p.plot(sds, np.log10(krumholz_ratio_model(sds, c=1, Z=0.5)), "r--",
       label="c=1, Z=0.5")
p.plot(sds, np.log10(krumholz_ratio_model(sds, c=4, Z=0.5)), "g-.",
       label="c=4, Z=0.5", linewidth=2)
p.plot(sds, np.log10(krumholz_ratio_model(sds, c=3, Z=0.5)), "m.",
       label="c=3, Z=0.5", linewidth=2)
p.plot(sds, np.log10(krumholz_ratio_model(sds, c=2, Z=1.0)), "b-",
       label="c=2, Z=1.0", linewidth=2)
p.legend(loc='lower right', frameon=True)
p.tight_layout()

save_name = "ratio_totalsigma_w_krumholzmodel_perpix_outliers"
p.savefig(proposal_figures_path("{}.pdf".format(save_name)))
p.savefig(proposal_figures_path("{}.png".format(save_name)))
p.close()

# Line Ratios

cube_reproj = spatial_reproject(cube, hi_cube.header)

co10_name = os.path.expanduser("~/MyRAID/M33/co/m33.combo.20.fits")
co10_cube = SpectralCube.read(co10_name)
co10_cube = co10_cube.with_spectral_unit(u.m / u.s,
                                         velocity_convention='radio',
                                         rest_value=115.271 * u.GHz)
# Spatial reprojection
co10_cube = spatial_reproject(co10_cube, hi_cube.header)
# noise_co10 = Noise(co10_cube)
# noise_co10.estimate_noise()
# noise_co10.get_scale_cube()
# masked_co10_cube = simple_masking(co10_cube, noise_co10.scale_cube, min_sig=2,
#                                   max_sig=5, min_pix=10)
# co10_mom0 = masked_co10_cube.moment0()

co32_name = os.path.expanduser("~/MyRAID/M33/co/m33.co32.fits")
co32_cube = SpectralCube.read(co32_name)
co32_cube = spatial_reproject(co32_cube, hi_cube.header)

# Manually assign the beam in
co32_cube.beam = Beam(18 * u.arcsec)
co32_cube.meta['beam'] = Beam(18 * u.arcsec)
co32_cube = co32_cube.convolve_to(co10_cube.beam)

co21_cube = cube.convolve_to(co10_cube.beam)

# noise_co32 = Noise(co32_cube)
# noise_co32.estimate_noise()
# noise_co32.get_scale_cube()
# masked_co32_cube = simple_masking(co32_cube, noise_co32.scale_cube, min_sig=2,
#                                   max_sig=5, min_pix=10)
# co32_mom0 = masked_co32_cube.moment0()


# Project everything to the HI grid. Then smooth to the CO(1-0)

# co10_mom0_reproj = co10_mom0.reproject(co_mom0.header)
# co32_mom0_reproj = co32_mom0.reproject(co_mom0.header)

# co32_mom0_reproj = co32_mom0_reproj.convolve_to(co10_cube.beam)

# Create stacked spectra for CO10, CO32

peak_co21_vels = \
    fits.open(iram_co21_data_path("m33.co21_iram.masked.peakvel.hireprojection.fits"))[0]
peak_co21_vels = Projection.from_hdu(peak_co21_vels)

co10_stacked = []
co21_stacked = []
co32_stacked = []

for num in [137, 139]:

    plane = (cloud_mask.data == num).sum(0) > 0

    # Now project the mask footprint onto the HI grid
    plane_reproj = reproject_interp((plane,
                                     WCS(cloud_mask.header).celestial),
                                    hi_cube.wcs.celestial,
                                    shape_out=hi_cube.shape[1:])[0] > 0

    xy_posns = np.where(plane_reproj)

    co10_shifted = []

    for y, x in zip(*xy_posns):
        co10_shifted.append(spectrum_shifter(co10_cube[:, y, x],
                                             peak_co21_vels[y, x], gal.vsys))

    stacked_co10 = np.zeros(co10_cube.shape[0])
    co10_count = np.zeros((co10_cube.shape[0]))

    for ii, spec in enumerate(co10_shifted):
        finites = np.isfinite(spec.value)
        stacked_co10[finites] += spec.value[finites]
        co10_count += finites

    average_co10 = stacked_co10 / co10_count
    co10_stacked.append(average_co10)

    co21_shifted = []

    for y, x in zip(*xy_posns):
        co21_shifted.append(spectrum_shifter(co21_cube[:, y, x],
                                             peak_co21_vels[y, x], gal.vsys))

    stacked_co21 = np.zeros(co21_cube.shape[0])
    co21_count = np.zeros((co21_cube.shape[0]))

    for ii, spec in enumerate(co21_shifted):
        finites = np.isfinite(spec.value)
        stacked_co21[finites] += spec.value[finites]
        co21_count += finites

    average_co21 = stacked_co21 / co21_count
    co21_stacked.append(average_co21)

    co32_shifted = []

    for y, x in zip(*xy_posns):
        co32_shifted.append(spectrum_shifter(co32_cube[:, y, x],
                                             peak_co21_vels[y, x], gal.vsys))

    stacked_co32 = np.zeros(co32_cube.shape[0])
    co32_count = np.zeros((co32_cube.shape[0]))

    for ii, spec in enumerate(co32_shifted):
        finites = np.isfinite(spec.value)
        stacked_co32[finites] += spec.value[finites]
        co32_count += finites

    average_co32 = stacked_co32 / co32_count
    co32_stacked.append(average_co32)

# Stacked spectra
fig, ax = p.subplots(2, 1,
                     sharex=True,
                     sharey=False, num=1)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)

fig.text(0.5, 0.02, 'Velocity (km/s)', ha='center')

ax[0].plot((hi_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
           hi_spectra[:, 136, 0].value, 'b',
           drawstyle='steps-mid', label='HI')
# For the legend
ax[0].plot([], [], 'k-', linewidth=3, label="CO(1-0)")
ax[0].plot([], [], 'r--', linewidth=3, label="CO(2-1)")
ax[0].plot([], [], 'g-.', linewidth=3, label="CO(3-2)")
ax2 = ax[0].twinx()
ax2.plot((cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
         co_spectra[:, 136, 0] * 1000.,
         # co21_shifted[0] * 1000.,
         'r--', linewidth=3,
         drawstyle='steps-mid', label='CO(2-1)')
ax2.plot((co10_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
         co10_stacked[0] * 1000.,
         'k-', linewidth=3,
         drawstyle='steps-mid', label='CO(1-0)', alpha=0.6)
ax2.plot((co32_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
         co32_stacked[0] * 1000.,
         'g-.', linewidth=3,
         drawstyle='steps-mid', label='CO(3-2)')
align_yaxis(ax[0], 0, ax2, 0)
# ax[0].set_ylim([0, 1])
ax[0].set_xlim([-40, 40])
ax[0].set_ylabel("HI Intensity (K)")
ax2.set_ylabel("CO Intensity (mK)")
ax[0].grid()
ax[0].text(-35, 30, "(Northern) Cloud 139", bbox=props)


ax[0].legend(frameon=True, loc='upper right')

ax[1].plot((hi_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
           hi_spectra[:, 138, 0].value, 'b',
           drawstyle='steps-mid', label='HI')
ax2 = ax[1].twinx()
ax2.plot((cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
         co_spectra[:, 138, 0] * 1000.,
         # co21_shifted[1] * 1000.,
         'r--', linewidth=3,
         drawstyle='steps-mid', label='CO(2-1)')
ax2.plot((co10_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
         co10_stacked[1] * 1000.,
         'k-', linewidth=3,
         drawstyle='steps-mid', label='CO(1-0)', alpha=0.6)
# ax2.plot((co32_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
#          co32_stacked[1] * 1000.,
#          'g-.',
#          drawstyle='steps-mid', label='CO(3-2)')
ax2.set_ylim(-20, 100)
align_yaxis(ax[1], 0, ax2, 0)
# ax[1].set_ylim([0, 1])
ax[1].set_xlim([-40, 40])
ax[1].set_ylabel("HI Intensity (K)")
ax2.set_ylabel("CO Intensity (mK)")
ax[1].grid()
ax[1].text(-35, 10, "(Southern) Cloud 137", bbox=props)

save_name = "clouds137_139_hi_co_lines"
p.savefig(proposal_figures_path("{}.pdf".format(save_name)))
p.savefig(proposal_figures_path("{}.png".format(save_name)))
p.close()

# Only CO(2-1)
fig, ax = p.subplots(2, 1,
                     sharex=True,
                     sharey=False, num=1)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)

fig.text(0.5, 0.02, 'Velocity (km/s)', ha='center')

ax[0].plot((hi_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
           hi_spectra[:, 136, 0].value, 'b',
           drawstyle='steps-mid', label='HI')
# For the legend
ax[0].plot([], [], 'r--', linewidth=3, label="CO(2-1)")
ax2 = ax[0].twinx()
ax2.plot((cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
         co_spectra[:, 136, 0] * 1000.,
         'r--', linewidth=3,
         drawstyle='steps-mid', label='CO(2-1)')
align_yaxis(ax[0], 0, ax2, 0)
# ax[0].set_ylim([0, 1])
ax[0].set_xlim([-40, 40])
ax[0].set_ylabel("HI Intensity (K)")
ax2.set_ylabel("CO Intensity (mK)")
ax[0].grid()
ax[0].text(-35, 30, "(Northern) Cloud 139", bbox=props)


ax[0].legend(frameon=True, loc='upper right')

ax[1].plot((hi_cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
           hi_spectra[:, 138, 0].value, 'b',
           drawstyle='steps-mid', label='HI')
ax2 = ax[1].twinx()
ax2.plot((cube.spectral_axis - gal.vsys).to(u.km / u.s).value,
         co_spectra[:, 138, 0] * 1000.,
         'r--', linewidth=3,
         drawstyle='steps-mid', label='CO(2-1)')
ax2.set_ylim(-20, 100)
align_yaxis(ax[1], 0, ax2, 0)
# ax[1].set_ylim([0, 1])
ax[1].set_xlim([-40, 40])
ax[1].set_ylabel("HI Intensity (K)")
ax2.set_ylabel("CO Intensity (mK)")
ax[1].grid()
ax[1].text(-35, 10, "(Southern) Cloud 137", bbox=props)

save_name = "clouds137_139_hi_co21"
p.savefig(proposal_figures_path("{}.pdf".format(save_name)))
p.savefig(proposal_figures_path("{}.png".format(save_name)))
p.close()


# Pick out the points belonging to each cloud

# def check_posns(xy_posns):
#     vals = []
#     for y, x in zip(*xy_posns):
#         posn = np.where(np.logical_and(tab["xpix"] == y, tab["ypix"] == x))
#         if posn[0].size != 0:
#             vals.append(posn[0][0])
#     return np.array(vals)


# num = 151

# plane = (cloud_mask.data == num).sum(0) > 0

# # Now project the mask footprint onto the HI grid
# plane_reproj = reproject_interp((plane,
#                                  WCS(cloud_mask.header).celestial),
#                                 hi_cube.wcs.celestial,
#                                 shape_out=hi_cube.shape[1:])[0] > 0

# xy_posns = np.where(plane_reproj)

# posns_151 = check_posns(xy_posns)
