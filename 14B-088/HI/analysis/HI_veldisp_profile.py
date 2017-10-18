import matplotlib.pyplot as p
from pandas import read_csv
from spectral_cube import SpectralCube, Projection
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle
from astropy.io import fits
import os
from os.path import join as osjoin

from cube_analysis.profiles import radial_profile

from paths import (fourteenB_HI_data_path, fourteenB_HI_file_dict,
                   fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, alltables_path)

from galaxy_params import gal_feath as gal
from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_twopanel_figure)

stack_figure_folder = allfigs_path("stacked_profiles")
if not os.path.exists(stack_figure_folder):
    os.mkdir(allfigs_path(stack_figure_folder))

props_figure_folder = allfigs_path("HI_properties")
if not os.path.exists(props_figure_folder):
    os.mkdir(allfigs_path(props_figure_folder))

lwidth_hdu = fits.open(fourteenB_HI_file_dict["LWidth"])[0]
lwidth = Projection.from_hdu(lwidth_hdu).to(u.km / u.s)

lwidth_feath_hdu = fits.open(fourteenB_wGBT_HI_file_dict["LWidth"])[0]
lwidth_feath = Projection.from_hdu(lwidth_feath_hdu).to(u.km / u.s)

# Create a radial profile of the HI vel disp out to 8 kpc.
# Beyond 8 kpc, noise is dominant. It may be some reflection of the
# warp, but I don't trust it right now.
rs, sd, sd_sigma = radial_profile(gal, lwidth, max_rad=8 * u.kpc)
rs_feath, sd_feath, sd_feath_sigma = \
    radial_profile(gal, lwidth_feath, max_rad=8 * u.kpc)
onecolumn_figure(font_scale=1.)
p.errorbar(rs.value, sd.value,
           yerr=sd_sigma.value, fmt="-",
           drawstyle='steps-mid', label='VLA')
p.errorbar(rs_feath.value, sd_feath.value,
           yerr=sd_feath_sigma.value, fmt="--",
           drawstyle='steps-mid', label='VLA + GBT')
p.xlabel("Radius (kpc)")
p.ylabel("HI Velocity Dispersion (km/s)")
p.legend(loc='lower left', frameon=True)
p.grid()
p.tight_layout()
p.savefig(osjoin(props_figure_folder, "hi_veldisp_profile.png"))
p.savefig(osjoin(props_figure_folder, "hi_veldisp_profile.pdf"))
p.close()

# Create the North and South portions.
rs_n, sd_n, sd_sigma_n = \
    radial_profile(gal, lwidth, max_rad=8 * u.kpc,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
rs_s, sd_s, sd_sigma_s = \
    radial_profile(gal, lwidth, max_rad=8 * u.kpc,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))
p.plot(rs_n.value, sd_n.value, "-.", label="North",
       drawstyle='steps-mid')
p.plot(rs_s.value, sd_s.value, "--", label="South",
       drawstyle='steps-mid')
p.errorbar(rs.value, sd.value,
           yerr=sd_sigma.value, fmt="-",
           drawstyle='steps-mid', label='Total')
p.xlabel("Radius (kpc)")
p.ylabel("HI Velocity Dispersion (km/s)")
p.grid()
p.legend(frameon=True)
p.tight_layout()
p.savefig(osjoin(props_figure_folder, "hi_veldisp_profile_n_s.png"))
p.savefig(osjoin(props_figure_folder, "hi_veldisp_profile_n_s.pdf"))
p.close()

rs_n, sd_n, sd_sigma_n = \
    radial_profile(gal, lwidth_feath, max_rad=8 * u.kpc,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
rs_s, sd_s, sd_sigma_s = \
    radial_profile(gal, lwidth_feath, max_rad=8 * u.kpc,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))
p.plot(rs_n.value, sd_n.value, "-.", label="North",
       drawstyle='steps-mid')
p.plot(rs_s.value, sd_s.value, "--", label="South",
       drawstyle='steps-mid')
p.errorbar(rs_feath.value, sd_feath.value,
           yerr=sd_feath_sigma.value, fmt="-",
           drawstyle='steps-mid', label='Total')
p.xlabel("Radius (kpc)")
p.ylabel("HI Velocity Dispersion (km/s)")
p.grid()
p.legend(frameon=True)
p.tight_layout()
p.savefig(osjoin(props_figure_folder, "hi_veldisp_profile_n_s_feather.png"))
p.savefig(osjoin(props_figure_folder, "hi_veldisp_profile_n_s_feather.pdf"))
p.close()


# Now load in the line stacking fits with the same bin size
hi_radial_fits = \
    read_csv(fourteenB_HI_data_path("tables/hi_gaussian_hwhm_totalprof_fits_radial.csv"))
twocolumn_twopanel_figure(font_scale=1.25)

fig, ax = p.subplots(1, 2, sharex=True, sharey=True)

ax[0].errorbar(rs.value, sd.value,
               yerr=sd_sigma.value, fmt="-",
               drawstyle='steps-mid', label='Avg. Line Width')
ax[0].set_xlabel("Radius (kpc)")
ax[0].set_ylabel("HI Velocity Dispersion (km/s)")
# Now the stacked fits
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['rotsub_sigma'],
               yerr=hi_radial_fits['rotsub_sigma_stderr'],
               fmt='D', label='Rot. Stack', alpha=0.5)
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['centsub_sigma'],
               yerr=hi_radial_fits['centsub_sigma_stderr'],
               fmt='o', label='Cent. Stack', alpha=0.5)
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_sigma'],
               yerr=hi_radial_fits['peaksub_sigma_stderr'],
               fmt='^', label='Peak Stack', alpha=0.5)
ax[0].grid()
ax[0].legend(frameon=True)
ax[0].text(0, 14, "VLA",
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[1].errorbar(rs.value, sd_feath.value,
               yerr=sd_feath_sigma.value, fmt="-",
               drawstyle='steps-mid', label='Averaged Line Width')
ax[1].set_xlabel("Radius (kpc)")
# Now the stacked fits
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['rotsub_sigma'],
               yerr=hi_radial_fits['rotsub_sigma_stderr'],
               fmt='D', label='Rot. Stack', alpha=0.5)
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['centsub_sigma'],
               yerr=hi_radial_fits['centsub_sigma_stderr'],
               fmt='o', label='Cent. Stack', alpha=0.5)
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_sigma'],
               yerr=hi_radial_fits['peaksub_sigma_stderr'],
               fmt='^', label='Peak Stack', alpha=0.5)
ax[1].grid()
ax[1].text(6.1, 14, "VLA + GBT",
           bbox={"boxstyle": "square", "facecolor": "w"})

p.tight_layout()

p.savefig(osjoin(stack_figure_folder, "hi_veldisp_w_stackedfits.png"))
p.savefig(osjoin(stack_figure_folder, "hi_veldisp_w_stackedfits.pdf"))
p.close()


# Let's compare the line width from the second moment to the Gaussian width
rot_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_100pc.fits"))
cent_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_100pc.fits"))
peakvel_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_100pc.fits"))
twocolumn_twopanel_figure(font_scale=1.2)
fig, ax = p.subplots(1, 3, sharey=True)
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['rotsub_sigma'],
               yerr=hi_radial_fits['rotsub_sigma_stderr'],
               fmt='D', label='Gaussian Fit', alpha=0.5)
ax[0].plot(hi_radial_fits['bin_center'],
           rot_stack.linewidth_sigma().value / 1000.,
           label="Moment")
ax[0].text(5, 11.5, "Rotation\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].legend(frameon=True, loc='lower right')
ax[0].grid()
# ax[0].set_xticklabels([])
ax[0].set_ylabel("HI Velocity Dispersion (km/s)")
ax[0].set_xlabel("Radius (kpc)")
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['centsub_sigma'],
               yerr=hi_radial_fits['centsub_sigma_stderr'],
               fmt='D', label='Gaussian Fit', alpha=0.5)
ax[1].plot(hi_radial_fits['bin_center'],
           cent_stack.linewidth_sigma().value / 1000.,
           label="Moment")
ax[1].text(5, 11.5, "Centroid\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[2].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_sigma'],
               yerr=hi_radial_fits['peaksub_sigma_stderr'],
               fmt='D', label='Gaussian Fit', alpha=0.5)
ax[2].plot(hi_radial_fits['bin_center'],
           cent_stack.linewidth_sigma().value / 1000.,
           label="Moment")
ax[2].text(5, 11.5, "Peak Vel.\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[2].grid()
ax[2].set_xlabel("Radius (kpc)")
p.tight_layout()
fig.savefig(osjoin(stack_figure_folder, "hi_veldisp_avg_vs_stackedfits.png"))
fig.savefig(osjoin(stack_figure_folder, "hi_veldisp_avg_vs_stackedfits.pdf"))
p.close()

rot_stack = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_100pc.fits"))
cent_stack = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_100pc.fits"))
peakvel_stack = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_100pc.fits"))
fig, ax = p.subplots(1, 3, sharey=True)
ax[0].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['rotsub_feath_sigma'],
               yerr=hi_radial_fits['rotsub_feath_sigma_stderr'],
               fmt='D', label='Gaussian Fit', alpha=0.5)
ax[0].plot(hi_radial_fits['bin_center'],
           rot_stack.linewidth_sigma().value / 1000.,
           label="Moment")
ax[0].text(4.75, 14.5, "Rotation\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[0].legend(frameon=True, loc='lower left')
ax[0].grid()
# ax[0].set_xticklabels([])
ax[0].set_ylabel("HI Velocity Dispersion (km/s)")
ax[0].set_xlabel("Radius (kpc)")
ax[1].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['centsub_feath_sigma'],
               yerr=hi_radial_fits['centsub_feath_sigma_stderr'],
               fmt='D', label='Gaussian Fit', alpha=0.5)
ax[1].plot(hi_radial_fits['bin_center'],
           cent_stack.linewidth_sigma().value / 1000.,
           label="Moment")
ax[1].text(4.75, 14.5, "Centroid\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[2].errorbar(hi_radial_fits['bin_center'],
               hi_radial_fits['peaksub_feath_sigma'],
               yerr=hi_radial_fits['peaksub_feath_sigma_stderr'],
               fmt='D', label='Gaussian Fit', alpha=0.5)
ax[2].plot(hi_radial_fits['bin_center'],
           cent_stack.linewidth_sigma().value / 1000.,
           label="Moment")
ax[2].text(4.75, 14.5, "Peak Vel.\nsubtracted",
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[2].grid()
ax[2].set_xlabel("Radius (kpc)")
ax[2].set_ylim([5, 17])
p.tight_layout()
fig.savefig(osjoin(stack_figure_folder, "hi_veldisp_avg_vs_stackedfits_feath.png"))
fig.savefig(osjoin(stack_figure_folder, "hi_veldisp_avg_vs_stackedfits_feath.pdf"))
p.close()

default_figure()
