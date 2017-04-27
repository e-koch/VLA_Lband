
'''
Exploring the properties of the skewness and kurtosis arrays.
'''

from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import numpy as np
import matplotlib.pyplot as p
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle

from analysis.paths import (fourteenB_HI_data_path, paper1_figures_path,
                            iram_co21_data_path, data_path)
from constants import (hi_freq, cube_name, moment0_name, lwidth_name,
                       skew_name, kurt_name, mask_name, moment1_name)
from plotting_styles import (twocolumn_figure, onecolumn_figure,
                             default_figure, twocolumn_twopanel_figure,
                             onecolumn_twopanel_figure,
                             onecolumn_Npanel_figure)

from HI_veldisp_profile import radial_profile
from HI_radial_profile import surfdens_radial_profile
from galaxy_params import gal


mom0_hdu = fits.open(fourteenB_HI_data_path(moment0_name))[0]
mom0 = Projection.from_hdu(mom0_hdu)

mom1_hdu = fits.open(fourteenB_HI_data_path(moment1_name))[0]
mom1 = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_HI_data_path(lwidth_name))[0]
lwidth = Projection.from_hdu(lwidth_hdu).to(u.km / u.s)

skew_hdu = fits.open(fourteenB_HI_data_path(skew_name))[0]
# Remove a few bad outliers from skew
skew_hdu.data[skew_hdu.data > 10] = np.NaN
skew_hdu.data[skew_hdu.data < -10] = np.NaN
skew = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_HI_data_path(kurt_name))[0]
kurt_hdu.data[kurt_hdu.data > 10] = np.NaN
kurt_hdu.data[kurt_hdu.data < -10] = np.NaN
kurt = Projection.from_hdu(kurt_hdu)


dr = 100 * u.pc

radii_skew, sdprof_skew, sdprof_sigma_skew = \
    radial_profile(gal, skew, dr=dr, max_rad=8 * u.kpc)

rs_skew_n, sd_skew_n, sd_skew_sigma_n = \
    radial_profile(gal, skew, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
rs_skew_s, sd_skew_s, sd_skew_sigma_s = \
    radial_profile(gal, skew, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))

radii_kurt, sdprof_kurt, sdprof_sigma_kurt = \
    radial_profile(gal, kurt, dr=dr, max_rad=8 * u.kpc)

rs_kurt_n, sd_kurt_n, sd_kurt_sigma_n = \
    radial_profile(gal, kurt, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
rs_kurt_s, sd_kurt_s, sd_kurt_sigma_s = \
    radial_profile(gal, kurt, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))

onecolumn_twopanel_figure(font_scale=1.2)

fig, ax = p.subplots(2, 1, sharex=True)
ax[0].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value,
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[0].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", )
ax[0].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", )
ax[0].set_ylim([-0.5, 0.4])
ax[0].legend(frameon=True, loc='lower right')
ax[0].grid()
# ax[0].set_xticklabels([])
ax[0].set_ylabel("Skewness")

ax[1].errorbar(radii_kurt.value, sdprof_kurt.value, yerr=sdprof_sigma_kurt.value,
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[1].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North")
ax[1].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South")
ax[1].set_ylim([-0.5, 0.4])

ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].set_ylabel("Kurtosis")

p.tight_layout()

fig.savefig(paper1_figures_path("hi_skew_kurt_profile_n_s.png"))
fig.savefig(paper1_figures_path("hi_skew_kurt_profile_n_s.pdf"))
p.close()

# Make a 4-panel stack of avg. profiles w/ mom0 and lwidth for comparison.

rs_lwidth, sd_lwidth, sd_lwidth_sigma = \
    radial_profile(gal, lwidth, max_rad=8 * u.kpc)

rs_lwidth_n, sd_lwidth_n, sd_lwidth_sigma_n = \
    radial_profile(gal, lwidth, max_rad=8 * u.kpc,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
rs_lwidth_s, sd_lwidth_s, sd_lwidth_sigma_s = \
    radial_profile(gal, lwidth, max_rad=8 * u.kpc,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))

rs_surf, sd_surf, sd_surf_sigma = \
    surfdens_radial_profile(gal, mom0=mom0, dr=dr,
                            restfreq=hi_freq, max_rad=8 * u.kpc)
rs_surf_n, sd_surf_n, sd_surf_sigma_n = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                            -0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq, max_rad=8 * u.kpc)
rs_surf_s, sd_surf_s, sd_surf_sigma_s = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq, max_rad=8 * u.kpc)

onecolumn_Npanel_figure(N=3, font_scale=1.0)

fig, ax = p.subplots(4, 1, sharex=True)

ax[0].errorbar(rs_surf.value, sd_surf.value,
               yerr=sd_surf_sigma.value,
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[0].errorbar(rs_surf_n.value, sd_surf_n.value, yerr=sd_surf_sigma_n.value,
               drawstyle='steps-mid', label="North", )
ax[0].errorbar(rs_surf_s.value, sd_surf_s.value, yerr=sd_surf_sigma_s.value,
               drawstyle='steps-mid', label="South", )
ax[0].set_ylim([2, 9])
ax[0].legend(frameon=True, loc='lower center')
ax[0].grid()
ax[0].set_ylabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")

ax[1].errorbar(rs_lwidth.value, sd_lwidth.value,
               yerr=sd_lwidth_sigma.value,
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[1].errorbar(rs_lwidth_n.value, sd_lwidth_n.value,
               yerr=sd_lwidth_sigma_n.value,
               drawstyle='steps-mid', label="North", )
ax[1].errorbar(rs_lwidth_s.value, sd_lwidth_s.value,
               yerr=sd_lwidth_sigma_s.value,
               drawstyle='steps-mid', label="South", )
ax[1].set_ylim([4, 14])
ax[1].grid()
ax[1].set_ylabel(r"Velocity Dispersion (km s$^{-1}$)")

ax[2].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value,
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[2].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", )
ax[2].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", )
ax[2].set_ylim([-0.5, 0.4])
ax[2].grid()
ax[2].set_ylabel("Skewness")

ax[3].errorbar(radii_kurt.value, sdprof_kurt.value,
               yerr=sdprof_sigma_kurt.value,
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[3].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North")
ax[3].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South")
ax[3].set_ylim([-0.5, 0.4])

ax[3].grid()
ax[3].set_xlabel("Radius (kpc)")
ax[3].set_ylabel("Kurtosis")

p.tight_layout()

fig.savefig(paper1_figures_path("hi_skew_kurt_sd_lwidth_profile_n_s.png"))
fig.savefig(paper1_figures_path("hi_skew_kurt_sd_lwidth_profile_n_s.pdf"))
p.close()

default_figure()
