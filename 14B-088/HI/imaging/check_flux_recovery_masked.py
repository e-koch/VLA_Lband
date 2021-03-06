
'''
Compare the flux in the cubes with their respective signal masks
'''

import numpy as np
from spectral_cube import SpectralCube
import os
import astropy.io.fits as fits
import astropy.units as u
import matplotlib.pyplot as plt
from pandas import DataFrame
import scipy.ndimage as nd

from cube_analysis.feather_cubes import flux_recovery

from paths import (fourteenB_HI_data_path, fourteenB_HI_file_dict,
                   fourteenB_wGBT_HI_file_dict,
                   gbt_HI_data_path, allfigs_path, paper1_tables_path)
from constants import pb_lim, hi_mass_conversion_Jy, distance
from plotting_styles import default_figure, onecolumn_figure

# Load the non-pb masked cube
vla_cube = SpectralCube.read(fourteenB_HI_file_dict["Cube"])
vla_source_mask = fits.open(fourteenB_HI_file_dict["Source_Mask"])[0].data > 0
vla_cube = vla_cube.with_mask(vla_source_mask)

gbt_name = gbt_HI_data_path("14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid_registered.fits")
gbt_cube = SpectralCube.read(gbt_name)

feathered_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["Cube"])
feather_source_mask = \
    fits.open(fourteenB_wGBT_HI_file_dict["Source_Mask"])[0].data > 0
feathered_cube = feathered_cube.with_mask(feather_source_mask)

pbcov = fits.open(fourteenB_HI_data_path("M33_14B-088_pbcov.fits"))[0]
mask = pbcov.data > pb_lim
# Cut to be the same shape as the PB masked cubes
mask = mask[nd.find_objects(mask)[0]]

total_vla_profile, total_gbt_profile = \
    flux_recovery(vla_cube, gbt_cube, mask=mask, num_cores=6)
total_feathered_profile, total_gbt_profile = \
    flux_recovery(feathered_cube, gbt_cube, mask=mask, num_cores=6)

vel_axis = vla_cube.spectral_axis.to(u.km / u.s).value

onecolumn_figure()

# Plot ratio b/w high-res to GBT total flux per channel
plt.plot(vel_axis, total_feathered_profile / total_gbt_profile,
         label='VLA + GBT')
plt.plot(vel_axis, total_vla_profile / total_gbt_profile, label="VLA",
         linestyle='--')
plt.ylim([0.15, 1.5])
plt.legend(frameon=True)
plt.grid(True)
plt.ylabel("VLA-to-GBT Flux Ratio")
plt.xlabel("Velocity (km / s)")
plt.tight_layout()
# raw_input("?")
plt.savefig(allfigs_path("Imaging/vla_gbt_masked_flux_recovery_ratio.png"))
plt.savefig(allfigs_path("Imaging/vla_gbt_masked_flux_recovery_ratio.pdf"))
plt.close()

# Plot the total spectra
plt.plot(vel_axis, total_gbt_profile,
         label='VLA + GBT')
plt.plot(vel_axis, total_vla_profile, label="Masked VLA",
         linestyle='--')
plt.plot(vel_axis, total_feathered_profile,
         label='Masked VLA + GBT', linestyle=":")
leg = plt.legend(frameon=True)
leg.get_frame().set_alpha(1.0)
plt.grid(True)
plt.ylim([-7, 80])
plt.xlim([-179 - 135, -179 + 137])
plt.ylabel("Total Flux (Jy)")
plt.xlabel("Velocity (km / s)")
bbox_kwarg = dict(facecolor='w', alpha=0.7)
plt.text(-258, 72, "North", ha='center', va='center', bbox=bbox_kwarg)
plt.text(-111, 72, "South", ha='center', va='center', bbox=bbox_kwarg)
plt.axvline(-179.6, color='k', linestyle='--', alpha=0.5, linewidth=3)
plt.fill_between([-73, -179 + 137], -7, 80, color='gray', alpha=0.5)
plt.text(-55, 37.5, "Milky Way HI", rotation='vertical',
         ha='center', va='center')
plt.tight_layout()

plt.savefig(allfigs_path("Imaging/vla_gbt_masked_flux_recovery.png"))
plt.savefig(allfigs_path("Imaging/vla_gbt_masked_flux_recovery.pdf"))
plt.close()

# We've summed up most of the data already. How about a mass estimate?
chan_width = np.abs(vel_axis[1] - vel_axis[0]) * u.km / u.s

vla_total_flux = np.sum(total_vla_profile) * chan_width
vla_mass = hi_mass_conversion_Jy * distance.to(u.Mpc)**2 * vla_total_flux

feathered_total_flux = np.sum(total_feathered_profile) * chan_width
feathered_mass = hi_mass_conversion_Jy * distance.to(u.Mpc)**2 * feathered_total_flux

gbt_total_flux = np.sum(total_gbt_profile) * chan_width
gbt_mass = hi_mass_conversion_Jy * distance.to(u.Mpc)**2 * gbt_total_flux


print("VLA HI Total Mass: {}".format(vla_mass))
print("GBT HI Total Mass: {}".format(gbt_mass))
print("VLA + GBT HI Total Mass: {}".format(feathered_mass))

df = DataFrame({"VLA Mass": [vla_mass.value],
                "GBT Mass": [gbt_mass.value],
                "VLA+GBT Mass": [feathered_mass.value]})
df.to_latex(paper1_tables_path("hi_masses_nomask.tex"))
df.to_csv(fourteenB_HI_data_wGBT_path("tables/hi_masses_withmask.csv",
                                      no_check=True))

# How much flux have we missed from avoiding the MW HI confusion regions?
# The last red-shifted side has ~6 Jy. Assume the shape will match the
# blue-shifted side. This occurs at -288.8 km/s.
miss_flux = total_gbt_profile[vel_axis < -288.8].sum() * chan_width

print("Missing flux from MW confusion: {}".format(miss_flux))

# Percentage of total flux?
miss_perc = miss_flux / (gbt_total_flux + miss_flux)
print("Percentage missing: {}".format(miss_perc))

default_figure()
