
import os
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

from paths import (fourteenB_HI_data_path, paper1_figures_path,
                   fourteenB_HI_file_dict, fourteenB_HI_data_wGBT_path,
                   fourteenB_wGBT_HI_file_dict, allfigs_path)
from plotting_styles import (default_figure, twocolumn_twopanel_figure,
                             twocolumn_figure)

'''
Plot the residual velocities from disk fit, and from subtracting the
smoothed version.
'''

direc = fourteenB_HI_data_path("diskfit_peakvels_noasymm_noradial_nowarp_output/",
                               no_check=True)

res = fits.open(os.path.join(direc, "rad.res.fits"))[0]
res_fit = fits.open(os.path.join(direc, "rad.fitres.fits"))[0]

default_figure()
ax2 = plt.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = plt.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

plt.savefig(allfigs_path("M33_residual_velocity_diskfit_fit.pdf"))
plt.savefig(allfigs_path("M33_residual_velocity_diskfit_fit.png"))

plt.close()

# Make a side-by-side centroid and residual plot.

twocolumn_twopanel_figure(fig_ratio=0.5)

ax2 = plt.subplot(122, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit.data / 1000., origin='lower', cmap='viridis',
                 vmax=20, vmin=-20)
ax2.set_xlabel("RA (J2000)")
ax2.coords[1].set_ticklabel_visible(False)

cb2 = plt.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

# moment1 = fits.open(fourteenB_HI_file_dict['Moment1'])[0]
peakvel = fits.open(fourteenB_HI_file_dict['PeakVels'])[0]

ax = plt.subplot(121, projection=WCS(peakvel.header))
im = ax.imshow(peakvel.data / 1000.,
               origin='lower', vmin=-300, vmax=-70,
               interpolation='nearest', cmap='viridis')

ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")

cbar = plt.colorbar(im)
cbar.set_label(r"Peak Velocity (km s$^{-1}$)")

plt.savefig(allfigs_path("M33_peakvel_residualfit.pdf"))
plt.savefig(allfigs_path("M33_peakvel_residualfit.png"))

plt.close()

# Now make similar plots with the feathered data

direc = fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_noradial_nowarp_output/",
                                    no_check=True)

res_feath = fits.open(os.path.join(direc, "rad.res.fits"))[0]
res_fit_feath = fits.open(os.path.join(direc, "rad.fitres.fits"))[0]

default_figure()
ax2 = plt.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit_feath.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = plt.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

plt.savefig(allfigs_path("M33_residual_velocity_diskfit_fit_feathered.pdf"))
plt.savefig(allfigs_path("M33_residual_velocity_diskfit_fit_feathered.png"))

plt.close()

# Make a side-by-side centroid and residual plot.

twocolumn_twopanel_figure(fig_ratio=0.5)

ax2 = plt.subplot(122, projection=WCS(res_fit_feath.header))
im2 = ax2.imshow(res_fit_feath.data / 1000., origin='lower', cmap='viridis',
                 vmax=20, vmin=-20)
ax2.set_xlabel("RA (J2000)")
ax2.coords[1].set_ticklabel_visible(False)

cb2 = plt.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

moment1 = fits.open(fourteenB_wGBT_HI_file_dict['Moment1'])[0]
peakvel = fits.open(fourteenB_wGBT_HI_file_dict['PeakVels'])[0]

ax = plt.subplot(121, projection=WCS(peakvel.header))
im = ax.imshow(peakvel.data / 1000.,
               origin='lower', vmin=-300, vmax=-70,
               interpolation='nearest', cmap='viridis')

ax.set_ylabel("DEC (J2000)")
ax.set_xlabel("RA (J2000)")

cbar = plt.colorbar(im)
cbar.set_label(r"Peak Velocity (km s$^{-1}$)")

plt.savefig(allfigs_path("M33_peakvel_residualfit_feathered.pdf"))
plt.savefig(allfigs_path("M33_peakvel_residualfit_feathered.png"))

plt.close()

default_figure()
