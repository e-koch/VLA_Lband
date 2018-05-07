
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
import os
from os.path import join as osjoin
import seaborn as sb
from astropy.table import Table, Column

from cube_analysis.profiles import (radial_profile, surfdens_radial_profile)

from paths import (fourteenB_HI_data_path, fourteenB_HI_file_dict,
                   fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path)
from constants import hi_freq, hi_mass_conversion
from plotting_styles import (twocolumn_figure, onecolumn_figure,
                             default_figure, twocolumn_twopanel_figure,
                             onecolumn_twopanel_figure,
                             onecolumn_Npanel_figure)

from galaxy_params import gal_feath as gal


prop_figure_folder = allfigs_path("HI_properties")
if not os.path.exists(prop_figure_folder):
    os.mkdir(allfigs_path(prop_figure_folder))

stack_figure_folder = allfigs_path("stacked_profiles")
if not os.path.exists(stack_figure_folder):
    os.mkdir(allfigs_path(stack_figure_folder))


# Color palette
cpal = sb.color_palette()

mom0_hdu = fits.open(fourteenB_HI_file_dict["Moment0"])[0]
mom0 = Projection.from_hdu(mom0_hdu)

mom1_hdu = fits.open(fourteenB_HI_file_dict["Moment1"])[0]
mom1 = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_HI_file_dict["LWidth"])[0]
lwidth = Projection.from_hdu(lwidth_hdu).to(u.km / u.s)

skew_hdu = fits.open(fourteenB_HI_file_dict["Skewness"])[0]
# Remove a few bad outliers from skew
skew_hdu.data[skew_hdu.data > 10] = np.NaN
skew_hdu.data[skew_hdu.data < -10] = np.NaN
skew = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_HI_file_dict["Kurtosis"])[0]
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

# Rename these so I can do a SD w/ and w/o comparison at the end
radii_skew_noSD = radii_skew.copy()
sdprof_skew_noSD = sdprof_skew.copy()
sdprof_sigma_skew_noSD = sdprof_sigma_skew.copy()

radii_skew_n_noSD = rs_skew_n.copy()
sdprof_skew_n_noSD = sd_skew_n.copy()
sdprof_sigma_skew_n_noSD = sd_skew_sigma_n.copy()

radii_skew_s_noSD = rs_skew_s.copy()
sdprof_skew_s_noSD = sd_skew_s.copy()
sdprof_sigma_skew_s_noSD = sd_skew_sigma_s.copy()

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
               yerr=sdprof_sigma_skew.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[0].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[0].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[0].set_ylim([-0.5, 0.4])
ax[0].legend(frameon=True, loc='lower right')
ax[0].grid()
ax[0].set_ylabel("Skewness")

ax[1].errorbar(radii_kurt.value, sdprof_kurt.value,
               yerr=sdprof_sigma_kurt.value,
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5, color='k',)
ax[1].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[1].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[1].set_ylim([-0.5, 0.4])

ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].set_ylabel("Kurtosis")

p.tight_layout()

fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_profile_n_s.png"))
fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_profile_n_s.pdf"))
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
                            restfreq=hi_freq, max_rad=8 * u.kpc,
                            mass_conversion=hi_mass_conversion)
rs_surf_n, sd_surf_n, sd_surf_sigma_n = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                            -0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq, max_rad=8 * u.kpc,
                            mass_conversion=hi_mass_conversion)
rs_surf_s, sd_surf_s, sd_surf_sigma_s = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq, max_rad=8 * u.kpc,
                            mass_conversion=hi_mass_conversion)

onecolumn_Npanel_figure(N=3, font_scale=1.0)

fig, ax = p.subplots(4, 1, sharex=True)

ax[0].errorbar(rs_surf.value, sd_surf.value,
               yerr=sd_surf_sigma.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[0].errorbar(rs_surf_n.value, sd_surf_n.value, yerr=sd_surf_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2])
ax[0].errorbar(rs_surf_s.value, sd_surf_s.value, yerr=sd_surf_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0])
ax[0].set_ylim([2, 9])
ax[0].legend(frameon=True, loc='lower center')
ax[0].grid()
ax[0].set_ylabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")

ax[1].errorbar(rs_lwidth.value, sd_lwidth.value,
               yerr=sd_lwidth_sigma.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[1].errorbar(rs_lwidth_n.value, sd_lwidth_n.value,
               yerr=sd_lwidth_sigma_n.value, color=cpal[2],
               drawstyle='steps-mid', label="North", )
ax[1].errorbar(rs_lwidth_s.value, sd_lwidth_s.value,
               yerr=sd_lwidth_sigma_s.value, color=cpal[0],
               drawstyle='steps-mid', label="South", )
ax[1].set_ylim([4, 14])
ax[1].grid()
ax[1].set_ylabel(r"Velocity Dispersion (km s$^{-1}$)")

ax[2].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[2].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2], )
ax[2].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0], )
ax[2].set_ylim([-0.5, 0.4])
ax[2].grid()
ax[2].set_ylabel("Skewness")

ax[3].errorbar(radii_kurt.value, sdprof_kurt.value,
               yerr=sdprof_sigma_kurt.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[3].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[3].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[3].set_ylim([-0.5, 0.4])

ax[3].grid()
ax[3].set_xlabel("Radius (kpc)")
ax[3].set_ylabel("Kurtosis")

p.tight_layout()

fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_sd_lwidth_profile_n_s.png"))
fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_sd_lwidth_profile_n_s.pdf"))
p.close()

# Now load in the stacked in the profiles and compare

cent_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_north_100pc.fits"))
cent_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_south_100pc.fits"))
cent_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_100pc.fits"))

peakvel_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_north_100pc.fits"))
peakvel_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_south_100pc.fits"))
peakvel_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_100pc.fits"))

rotation_stack_n = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_north_100pc.fits"))
rotation_stack_s = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_south_100pc.fits"))
rotation_stack = SpectralCube.read(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_100pc.fits"))


def skewness(cube):
    return cube.moment(order=3) / cube.linewidth_sigma()**3


def kurtosis(cube):
    return (cube.moment(order=4) / cube.linewidth_sigma()**4) - 3


dr = 100 * u.pc
max_radius = (8.0 * u.kpc).to(u.pc)

nbins = np.int(np.floor(max_radius / dr))
inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)
bin_centers = 0.5 * (inneredge + outeredge).to(u.kpc)

onecolumn_twopanel_figure(font_scale=1.2)

fig, ax = p.subplots(2, 1, sharex=True)
ax[0].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value, color=cpal[1],
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[0].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2], )
ax[0].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0], )
ax[0].plot(bin_centers, skewness(cent_stack).value, label='Total Stacked')
ax[0].plot(bin_centers, skewness(cent_stack_n).value, label='North Stacked')
ax[0].plot(bin_centers, skewness(cent_stack_s).value, label='South Stacked')
# ax[0].set_ylim([-0.5, 0.4])
ax[0].legend(frameon=True, loc='lower right')
ax[0].grid()
# ax[0].set_xticklabels([])
ax[0].set_ylabel("Skewness")

ax[1].errorbar(radii_kurt.value, sdprof_kurt.value,
               yerr=sdprof_sigma_kurt.value, color=cpal[1],
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[1].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[1].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[1].plot(bin_centers, kurtosis(cent_stack).value, label='Total Stacked')
ax[1].plot(bin_centers, kurtosis(cent_stack_n).value, label='North Stacked')
ax[1].plot(bin_centers, kurtosis(cent_stack_s).value, label='South Stacked')
# ax[1].set_ylim([-0.5, 0.4])

ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].set_ylabel("Kurtosis")

p.tight_layout()

fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_profile_n_s_w_cent_stacked.png"))
fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_profile_n_s_w_cent_stacked.pdf"))
p.close()

onecolumn_twopanel_figure(font_scale=1.0)

fig, ax = p.subplots(2, 1, sharex=True)
ax[0].plot(bin_centers, skewness(cent_stack).value, label='Total Cent.')
ax[0].plot(bin_centers, skewness(cent_stack_n).value, label='North Cent.')
ax[0].plot(bin_centers, skewness(cent_stack_s).value, label='South Cent.')

ax[0].plot(bin_centers, skewness(cent_stack).value, label='Total Peak',
           linestyle='--')
ax[0].plot(bin_centers, skewness(cent_stack_n).value, label='North Peak',
           linestyle='--')
ax[0].plot(bin_centers, skewness(cent_stack_s).value, label='South Peak',
           linestyle='--')

ax[0].plot(bin_centers, skewness(rotation_stack).value, label='Total Rot.',
           linestyle='-.')
ax[0].plot(bin_centers, skewness(rotation_stack_n).value, label='North Rot.',
           linestyle='-.')
ax[0].plot(bin_centers, skewness(rotation_stack_s).value, label='South Rot.',
           linestyle='-.')
ax[0].set_xlim([0, 12])
# ax[0].set_ylim([-0.5, 0.4])
ax[0].legend(frameon=True, loc='lower right')
ax[0].grid()
# ax[0].set_xticklabels([])
ax[0].set_ylabel("Skewness")

ax[1].plot(bin_centers, kurtosis(cent_stack).value, label='Total Cent.')
ax[1].plot(bin_centers, kurtosis(cent_stack_n).value, label='North Cent.')
ax[1].plot(bin_centers, kurtosis(cent_stack_s).value, label='South Cent.')

ax[1].plot(bin_centers, kurtosis(cent_stack).value, label='Total Peak',
           linestyle='--')
ax[1].plot(bin_centers, kurtosis(cent_stack_n).value, label='North Peak',
           linestyle='--')
ax[1].plot(bin_centers, kurtosis(cent_stack_s).value, label='South Peak',
           linestyle='--')

ax[1].plot(bin_centers, kurtosis(rotation_stack).value, label='Total Rot.',
           linestyle='-.')
ax[1].plot(bin_centers, kurtosis(rotation_stack_n).value, label='North Rot.',
           linestyle='-.')
ax[1].plot(bin_centers, kurtosis(rotation_stack_s).value, label='South Rot.',
           linestyle='-.')

ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].set_ylabel("Kurtosis")

p.tight_layout()

fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_stacked_comparison.png"))
fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_stacked_comparison.pdf"))
p.close()

# Now the feathered version
mom0_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]
mom0 = Projection.from_hdu(mom0_hdu)

# Beam handling still not finished. Remove per beam if it appears. The
# units will be wrong in the intermediate steps, but will be converted
# correctly at the end.
if "beam" in mom0.unit.to_string():
    mom0 = mom0 * u.beam

mom1_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Moment1"])[0]
mom1 = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_wGBT_HI_file_dict["LWidth"])[0]
lwidth = Projection.from_hdu(lwidth_hdu).to(u.km / u.s)

skew_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Skewness"])[0]
# Remove a few bad outliers from skew
skew_hdu.data[skew_hdu.data > 10] = np.NaN
skew_hdu.data[skew_hdu.data < -10] = np.NaN
skew = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Kurtosis"])[0]
kurt_hdu.data[kurt_hdu.data > 10] = np.NaN
kurt_hdu.data[kurt_hdu.data < -10] = np.NaN
kurt = Projection.from_hdu(kurt_hdu)

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
                            restfreq=hi_freq, max_rad=8 * u.kpc,
                            mass_conversion=hi_mass_conversion)
rs_surf_n, sd_surf_n, sd_surf_sigma_n = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                            -0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq, max_rad=8 * u.kpc,
                            mass_conversion=hi_mass_conversion)
rs_surf_s, sd_surf_s, sd_surf_sigma_s = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq, max_rad=8 * u.kpc,
                            mass_conversion=hi_mass_conversion)

# Save the feathered radial profiles as a table
tab_out = Table([Column(rs_surf, name="rad_bin"),
                 Column(sd_surf, name='Sigma_HI'),
                 Column(sd_surf_sigma, name='Sigma_HI_std'),
                 Column(sd_surf_n, name='Sigma_HI_north'),
                 Column(sd_surf_sigma_n, name='Sigma_HI_std_north'),
                 Column(sd_surf_s, name='Sigma_HI_south'),
                 Column(sd_surf_sigma_s, name='Sigma_HI_std_south'),
                 Column(sd_lwidth, name='lwidth'),
                 Column(sd_lwidth_sigma, name='lwidth_std'),
                 Column(sd_lwidth_n, name='lwidth_north'),
                 Column(sd_lwidth_sigma_n, name='lwidth_std_north'),
                 Column(sd_lwidth_s, name='lwidth_south'),
                 Column(sd_lwidth_sigma_s, name='lwidth_std_south'),
                 Column(sdprof_kurt, name='kurtosis'),
                 Column(sdprof_sigma_kurt, name='kurtosis_std'),
                 Column(sd_kurt_n, name='kurtosis_north'),
                 Column(sd_kurt_sigma_n, name='kurtosis_std_north'),
                 Column(sd_kurt_s, name='kurtosis_south'),
                 Column(sd_kurt_sigma_s, name='kurtosis_std_south'),
                 Column(sdprof_skew, name='skewness'),
                 Column(sdprof_sigma_skew, name='skewness_std'),
                 Column(sd_skew_n, name='skewness_north'),
                 Column(sd_skew_sigma_n, name='skewness_std_north'),
                 Column(sd_skew_s, name='skewness_south'),
                 Column(sd_skew_sigma_s, name='skewness_std_south'),
                 ])
tab_out.write(fourteenB_HI_data_wGBT_path("tables/moments_radial_profile_{0}{1}.fits"
                                          .format(int(dr.value), dr.unit),
                                          no_check=True),
              overwrite=True)

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
ax[0].set_ylim([3.5, 11.2])
ax[0].legend(frameon=True, loc='lower center')
ax[0].grid()
ax[0].set_ylabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")

ax[1].errorbar(rs_lwidth.value, sd_lwidth.value,
               yerr=sd_lwidth_sigma.value, color=cpal[1],
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[1].errorbar(rs_lwidth_n.value, sd_lwidth_n.value,
               yerr=sd_lwidth_sigma_n.value, color=cpal[2],
               drawstyle='steps-mid', label="North", )
ax[1].errorbar(rs_lwidth_s.value, sd_lwidth_s.value,
               yerr=sd_lwidth_sigma_s.value, color=cpal[0],
               drawstyle='steps-mid', label="South", )
ax[1].set_ylim([8, 16])
ax[1].grid()
ax[1].set_ylabel(r"Velocity Dispersion (km s$^{-1}$)")

ax[2].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value, color=cpal[1],
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[2].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[2].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[2].set_ylim([-0.7, 0.7])
ax[2].grid()
ax[2].set_ylabel("Skewness")

ax[3].errorbar(radii_kurt.value, sdprof_kurt.value,
               yerr=sdprof_sigma_kurt.value, color=cpal[1],
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[3].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[3].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[3].set_ylim([-0.6, 0.75])
ax[3].grid()
ax[3].set_xlabel("Radius (kpc)")
ax[3].set_ylabel("Kurtosis")

ax[0].set_xlim([0, 8])

p.tight_layout()

fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_sd_lwidth_profile_n_s_feather.png"))
fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_sd_lwidth_profile_n_s_feather.pdf"))
p.close()

# Feathered skew/kurt profiles in the N and S

onecolumn_Npanel_figure(N=2, font_scale=1.2)

fig, ax = p.subplots(2, 1, sharex=True)

ax[0].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=1.)
ax[0].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],
               linestyle='--')
ax[0].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],
               linestyle='-.')
ax[0].fill_between([0., 0.5], -0.75, 0.8, color='gray', alpha=0.5)
ax[0].set_ylim([-0.75, 0.7])
ax[0].grid()
ax[0].legend(frameon=True)
ax[0].axhline(0., color='gray', linestyle='dashed', linewidth=3,
              alpha=0.7, zorder=-2)
ax[0].set_ylabel("Skewness")

ax[1].errorbar(radii_kurt.value, sdprof_kurt.value,
               yerr=sdprof_sigma_kurt.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=1.)
ax[1].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],
               linestyle='--')
ax[1].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],
               linestyle='-.')
ax[1].fill_between([0., 0.5], -0.75, 0.8, color='gray', alpha=0.5)
ax[1].axhline(0., color='gray', linestyle='dashed', linewidth=3,
              alpha=0.7, zorder=-2)
ax[1].set_ylim([-0.6, 0.75])
ax[1].set_xlim([0, 8])
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].set_ylabel("Excess Kurtosis")

p.tight_layout()

fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_profile_n_s_feather.png"))
fig.savefig(osjoin(prop_figure_folder, "hi_skew_kurt_profile_n_s_feather.pdf"))
p.close()

cent_stack_n = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_north_100pc.fits"))
cent_stack_s = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_south_100pc.fits"))
cent_stack = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_100pc.fits"))

peakvel_stack_n = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_north_100pc.fits"))
peakvel_stack_s = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_south_100pc.fits"))
peakvel_stack = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_100pc.fits"))

rotation_stack_n = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_north_100pc.fits"))
rotation_stack_s = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_south_100pc.fits"))
rotation_stack = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_100pc.fits"))

onecolumn_twopanel_figure(font_scale=1.2)

fig, ax = p.subplots(2, 1, sharex=True)
ax[0].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value, color=cpal[1],
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[0].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[0].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[0].plot(bin_centers, skewness(cent_stack).value, label='Total Stacked')
ax[0].plot(bin_centers, skewness(cent_stack_n).value, label='North Stacked')
ax[0].plot(bin_centers, skewness(cent_stack_s).value, label='South Stacked')
ax[0].legend(frameon=True, loc='upper right')
ax[0].grid()
ax[0].set_ylabel("Skewness")

ax[1].errorbar(radii_kurt.value, sdprof_kurt.value,
               yerr=sdprof_sigma_kurt.value, color=cpal[1],
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[1].errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
               drawstyle='steps-mid', label="North", color=cpal[2],)
ax[1].errorbar(rs_kurt_s.value, sd_kurt_s.value, yerr=sd_kurt_sigma_s.value,
               drawstyle='steps-mid', label="South", color=cpal[0],)
ax[1].plot(bin_centers, kurtosis(cent_stack).value, label='Total Stacked')
ax[1].plot(bin_centers, kurtosis(cent_stack_n).value, label='North Stacked')
ax[1].plot(bin_centers, kurtosis(cent_stack_s).value, label='South Stacked')
ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].set_ylabel("Kurtosis")

ax[0].set_xlim([0, 14])

p.tight_layout()
fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_profile_n_s_w_cent_stacked_feather.png"))
fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_profile_n_s_w_cent_stacked_feather.pdf"))
p.close()

onecolumn_twopanel_figure(font_scale=1.0)

fig, ax = p.subplots(2, 1, sharex=True)
ax[0].plot(bin_centers, skewness(cent_stack).value, label='Total Cent.')
ax[0].plot(bin_centers, skewness(cent_stack_n).value, label='North Cent.')
ax[0].plot(bin_centers, skewness(cent_stack_s).value, label='South Cent.')

ax[0].plot(bin_centers, skewness(cent_stack).value, label='Total Peak',
           linestyle='--')
ax[0].plot(bin_centers, skewness(cent_stack_n).value, label='North Peak',
           linestyle='--')
ax[0].plot(bin_centers, skewness(cent_stack_s).value, label='South Peak',
           linestyle='--')

ax[0].plot(bin_centers, skewness(rotation_stack).value, label='Total Rot.',
           linestyle='-.')
ax[0].plot(bin_centers, skewness(rotation_stack_n).value, label='North Rot.',
           linestyle='-.')
ax[0].plot(bin_centers, skewness(rotation_stack_s).value, label='South Rot.',
           linestyle='-.')
ax[0].set_xlim([0, 12])
# ax[0].set_ylim([-0.5, 0.4])
ax[0].legend(frameon=True, loc='lower right')
ax[0].grid()
# ax[0].set_xticklabels([])
ax[0].set_ylabel("Skewness")

ax[1].plot(bin_centers, kurtosis(cent_stack).value, label='Total Cent.')
ax[1].plot(bin_centers, kurtosis(cent_stack_n).value, label='North Cent.')
ax[1].plot(bin_centers, kurtosis(cent_stack_s).value, label='South Cent.')

ax[1].plot(bin_centers, kurtosis(cent_stack).value, label='Total Peak',
           linestyle='--')
ax[1].plot(bin_centers, kurtosis(cent_stack_n).value, label='North Peak',
           linestyle='--')
ax[1].plot(bin_centers, kurtosis(cent_stack_s).value, label='South Peak',
           linestyle='--')

ax[1].plot(bin_centers, kurtosis(rotation_stack).value, label='Total Rot.',
           linestyle='-.')
ax[1].plot(bin_centers, kurtosis(rotation_stack_n).value, label='North Rot.',
           linestyle='-.')
ax[1].plot(bin_centers, kurtosis(rotation_stack_s).value, label='South Rot.',
           linestyle='-.')

ax[1].grid()
ax[1].set_xlabel("Radius (kpc)")
ax[1].set_ylabel("Kurtosis")

p.tight_layout()

fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_stacked_comparison_feather.png"))
fig.savefig(osjoin(stack_figure_folder, "hi_skew_kurt_stacked_comparison_feather.pdf"))
p.close()


# Finally, make a 2-panel plot of just the skewness profiles to highlight the
# N/S features
twocolumn_twopanel_figure()

fig, ax = p.subplots(1, 2, sharex=True, sharey=True)

ax[0].errorbar(radii_skew_noSD.value, sdprof_skew_noSD.value,
               yerr=sdprof_sigma_skew_noSD.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[0].errorbar(radii_skew_n_noSD.value, sdprof_skew_n_noSD.value,
               yerr=sdprof_sigma_skew_n_noSD.value, color=cpal[2],
               drawstyle='steps-mid', label="North",
               linestyle='--')
ax[0].errorbar(radii_skew_s_noSD.value, sdprof_skew_s_noSD.value,
               yerr=sdprof_sigma_skew_s_noSD.value, color=cpal[0],
               drawstyle='steps-mid', label="South",
               linestyle='-.')
ax[0].fill_between([0., 0.5], -0.75, 0.8, color='gray', alpha=0.5)
ax[0].set_ylim([-0.7, 0.78])
ax[0].set_xlim([0, 8])
ax[0].grid()
ax[0].legend(frameon=True, loc='upper center')
ax[0].text(7.5, 0.65, "VLA",
           horizontalalignment='right',
           verticalalignment='top',
           bbox={"boxstyle": "square", "facecolor": "w"})

ax[0].set_ylabel("Skewness")
ax[0].set_xlabel("Radius (kpc)")

ax[1].errorbar(radii_skew.value, sdprof_skew.value,
               yerr=sdprof_sigma_skew.value, color='k',
               drawstyle='steps-mid', label="Total", zorder=-1,
               alpha=0.5)
ax[1].errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
               drawstyle='steps-mid', label="North",
               linestyle='--', color=cpal[2],)
ax[1].errorbar(rs_skew_s.value, sd_skew_s.value, yerr=sd_skew_sigma_s.value,
               drawstyle='steps-mid', label="South",
               linestyle='-.', color=cpal[0],)
ax[1].fill_between([0., 0.5], -0.75, 0.8, color='gray', alpha=0.5)
ax[1].set_ylim([-0.7, 0.78])
ax[1].set_xlim([0, 8])
ax[1].grid()
ax[1].text(7.5, 0.65, "VLA+GBT",
           horizontalalignment='right',
           verticalalignment='top',
           bbox={"boxstyle": "square", "facecolor": "w"})
ax[1].set_xlabel("Radius (kpc)")

p.tight_layout()

fig.savefig(osjoin(prop_figure_folder, "hi_skew_profile_both_n_s.png"))
fig.savefig(osjoin(prop_figure_folder, "hi_skew_profile_both_n_s.pdf"))
p.close()

default_figure()
