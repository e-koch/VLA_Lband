
'''
Script to compare the flux recovered in the interferometer cube, the SD cube,
and (optionally) the feathered cube.
'''

import numpy as np
from spectral_cube import SpectralCube
import os
import astropy.io.fits as fits
import astropy.units as u
import matplotlib.pyplot as plt
from pandas import DataFrame

from cube_analysis.feather_cubes import flux_recovery

from paths import (seventeenB_HI_data_02kms_path,
                   seventeenB_HI_data_02kms_wGBT_path,
                   seventeenB_HI_data_1kms_path,
                   seventeenB_HI_data_1kms_wGBT_path,
                   data_path, allfigs_path, paper1_tables_path)
from constants import pb_lim, hi_mass_conversion_Jy, distance
from plotting_styles import default_figure, onecolumn_figure

# Set which of the cubes to feather
run_gbt_02kms = True
run_gbt_1kms = False

num_cores = 4
chunk = 8

if run_gbt_02kms:

    # Load the non-pb masked cube
    vla_cube = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.fits"))

    gbt_path = os.path.join(data_path, "GBT")
    gbt_name = os.path.join(gbt_path, "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_02kms_spectralregrid.fits")
    gbt_cube = SpectralCube.read(gbt_name)

    feathered_cube = SpectralCube.read(seventeenB_HI_data_02kms_wGBT_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.GBT_feathered.fits"))

    pbcov = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.pb.fits"))
    mask = pbcov[0].value > pb_lim

    total_vla_profile, total_gbt_profile = \
        flux_recovery(vla_cube, gbt_cube, mask=mask, num_cores=num_cores,
                      chunk=chunk)
    total_feathered_profile, total_gbt_profile = \
        flux_recovery(feathered_cube, gbt_cube, mask=mask, num_cores=num_cores,
                      chunk=chunk)

    vel_axis = vla_cube.spectral_axis.to(u.km / u.s).value

    onecolumn_figure()

    # Plot ratio b/w high-res to GBT total flux per channel
    plt.plot(vel_axis, total_feathered_profile / total_gbt_profile,
             label='VLA + GBT')
    plt.plot(vel_axis, total_vla_profile / total_gbt_profile, label="VLA",
             linestyle='--')
    # plt.axhline(1, zorder=-1, linestyle='--', color='b', alpha=0.5)
    plt.ylim([0.15, 1.5])
    plt.legend(frameon=True)
    plt.grid(True)
    plt.ylabel("VLA-to-GBT Flux Ratio")
    plt.xlabel("Velocity (km / s)")
    plt.tight_layout()
    plt.savefig(allfigs_path("Imaging/vla_gbt_17B_flux_recovery_ratio.png"))
    plt.savefig(allfigs_path("Imaging/vla_gbt_17B_flux_recovery_ratio.pdf"))
    plt.close()

    # Plot the total spectra
    plt.plot(vel_axis, total_gbt_profile,
             label='GBT')
    plt.plot(vel_axis, total_vla_profile, label="VLA",
             linestyle='--')
    plt.plot(vel_axis, total_feathered_profile,
             label='VLA + GBT', linestyle=":")
    plt.legend(frameon=True)
    plt.grid(True)
    plt.ylim([-3, 70])
    plt.ylabel("Total Flux (Jy)")
    plt.xlabel("Velocity (km / s)")
    plt.tight_layout()
    plt.savefig(allfigs_path("Imaging/vla_gbt_17B_flux_recovery.png"))
    plt.savefig(allfigs_path("Imaging/vla_gbt_17B_flux_recovery.pdf"))
    plt.close()

    # We've summed up most of the data already. How about a mass estimate?
    chan_width = np.abs(vel_axis[1] - vel_axis[0]) * u.km / u.s

    vla_total_flux = np.sum(total_vla_profile) * chan_width
    vla_mass = hi_mass_conversion_Jy * distance.to(u.Mpc)**2 * vla_total_flux

    feathered_total_flux = np.sum(total_feathered_profile) * chan_width
    feathered_mass = hi_mass_conversion_Jy * distance.to(u.Mpc)**2 * \
        feathered_total_flux

    gbt_total_flux = np.sum(total_gbt_profile) * chan_width
    gbt_mass = hi_mass_conversion_Jy * distance.to(u.Mpc)**2 * gbt_total_flux

    print("VLA HI Total Mass: {}".format(vla_mass))
    print("GBT HI Total Mass: {}".format(gbt_mass))
    print("VLA + GBT HI Total Mass: {}".format(feathered_mass))

    df = DataFrame({"VLA Mass": [vla_mass.value],
                    "GBT Mass": [gbt_mass.value],
                    "VLA+GBT Mass": [feathered_mass.value]})
    df.to_csv(seventeenB_HI_data_02kms_wGBT_path("tables/hi_masses_nomask.csv",
                                                 no_check=True))

default_figure()
