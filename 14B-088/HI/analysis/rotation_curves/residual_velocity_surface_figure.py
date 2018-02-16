
from os.path import join as osjoin
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from aplpy import FITSFigure

from paths import (fourteenB_HI_data_path, paper1_figures_path,
                   fourteenB_HI_file_dict, fourteenB_HI_data_wGBT_path,
                   fourteenB_wGBT_HI_file_dict, allfigs_path)
from plotting_styles import (default_figure, twocolumn_twopanel_figure,
                             twocolumn_figure)

'''
Plot the residual velocities from disk fit, and from subtracting the
smoothed version.
'''

fig_folder = allfigs_path("HI_maps")

direc = fourteenB_HI_data_path("diskfit_peakvels_noasymm_noradial_nowarp_output/",
                               no_check=True)

res = fits.open(osjoin(direc, "rad.res.fits"))[0]
res_fit = fits.open(osjoin(direc, "rad.fitres.fits"))[0]

default_figure()
ax2 = plt.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = plt.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

plt.savefig(osjoin(fig_folder, "M33_residual_velocity_diskfit_fit.pdf"))
plt.savefig(osjoin(fig_folder, "M33_residual_velocity_diskfit_fit.png"))

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

plt.savefig(osjoin(fig_folder, "M33_peakvel_residualfit.pdf"))
plt.savefig(osjoin(fig_folder, "M33_peakvel_residualfit.png"))

plt.close()

# # Now make similar plots with the feathered data

direc = fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_noradial_nowarp_output/",
                                    no_check=True)

res_feath = fits.open(osjoin(direc, "rad.res.fits"))[0]
res_fit_feath = fits.open(osjoin(direc, "rad.fitres.fits"))[0]

default_figure()
ax2 = plt.subplot(111, projection=WCS(res_fit.header))
im2 = ax2.imshow(res_fit_feath.data / 1000., origin='lower', cmap='seismic',
                 vmax=30, vmin=-30)

cb2 = plt.colorbar(im2)
cb2.set_label("Residual Velocity (km s$^{-1}$)")

plt.savefig(osjoin(fig_folder, "M33_residual_velocity_diskfit_fit_feathered.pdf"))
plt.savefig(osjoin(fig_folder, "M33_residual_velocity_diskfit_fit_feathered.png"))

plt.close()

# Make a side-by-side centroid and residual plot. And the line width map!

moment1 = fits.open(fourteenB_wGBT_HI_file_dict['Moment1'])[0]
peakvel = fits.open(fourteenB_wGBT_HI_file_dict['PeakVels'])[0]
lwidth = fits.open(fourteenB_wGBT_HI_file_dict['LWidth'])[0]

twocolumn_twopanel_figure(fig_ratio=0.5)

fig = plt.figure()

fig2 = FITSFigure(fits.PrimaryHDU(res_fit_feath.data / 1000.,
                                  header=res_fit_feath.header), figure=fig,
                  subplot=(1, 3, 2))
fig2.show_colorscale(cmap='coolwarm', vmax=13, vmin=-13)
fig2.hide_axis_labels()
fig2.hide_ytick_labels()
fig2.tick_labels.set_xformat('hh:mm')

fig2.add_colorbar(location='top')
fig2.colorbar.set_axis_label_text('Residual Velocity (km / s)')
fig2.colorbar.set_axis_label_font(size=10)
fig2.colorbar.set_font(size=10)
fig2.colorbar.set_ticks([-20, -10, 0, 10, 20])

fig1 = FITSFigure(fits.PrimaryHDU(peakvel.data / 1000.,
                                  header=peakvel.header), figure=fig,
                  subplot=(1, 3, 1))
fig1.show_colorscale(cmap='viridis', vmin=-300, vmax=-70)
fig1.hide_axis_labels()
# fig1.hide_ytick_labels()

fig1.add_colorbar(location='top')
fig1.colorbar.set_axis_label_text('Peak Velocity (km / s)')
fig1.colorbar.set_axis_label_font(size=10)
fig1.colorbar.set_font(size=10)
fig1.tick_labels.set_xformat('hh:mm')
fig1.tick_labels.set_yformat('dd:mm')

fig3 = FITSFigure(fits.PrimaryHDU(lwidth.data / 1000.,
                                  header=peakvel.header), figure=fig,
                  subplot=(1, 3, 3))
fig3.show_colorscale(cmap='viridis', vmin=5, vmax=20)
fig3.hide_axis_labels()
fig3.hide_ytick_labels()

fig3.add_colorbar(location='top')
fig3.colorbar.set_axis_label_text('Line Width (km / s)')
fig3.colorbar.set_axis_label_font(size=10)
fig3.colorbar.set_font(size=10)
fig3.tick_labels.set_xformat('hh:mm')

plt.tight_layout()

plt.savefig(osjoin(fig_folder, "M33_peakvel_residualfit_feathered.pdf"))
plt.savefig(osjoin(fig_folder, "M33_peakvel_residualfit_feathered.png"))

plt.close()

default_figure()
