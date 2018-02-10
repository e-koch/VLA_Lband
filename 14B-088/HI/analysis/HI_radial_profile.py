
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.coordinates import Angle
import matplotlib.pyplot as p
import matplotlib.patches as patches
from astropy.table import Table
from spectral_cube.lower_dimensional_structures import Projection
from astropy.io import fits
import seaborn as sb

import os

from cube_analysis.profiles import surfdens_radial_profile

from paths import (arecibo_HI_data_path, gbt_HI_data_path,
                   c_hi_analysispath, allfigs_path,
                   data_path, fourteenB_wGBT_HI_file_dict,
                   fourteenB_HI_file_dict)

from constants import hi_freq, hi_mass_conversion
from galaxy_params import gal
from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_twopanel_figure)

mom0_hdu = fits.open(fourteenB_HI_file_dict["Moment0"])[0]
mom0 = Projection.from_hdu(mom0_hdu)

# And the feathered version
mom0_feath_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]
mom0_feath = Projection.from_hdu(mom0_feath_hdu)

if not os.path.exists(allfigs_path("HI_properties")):
    os.mkdir(allfigs_path("HI_properties"))

# Bin size in pc
dr = 100 * u.pc

# Create a radial profile of HI
rs, sd, sd_sigma = surfdens_radial_profile(gal, mom0=mom0, dr=dr,
                                           restfreq=hi_freq,
                                           mass_conversion=hi_mass_conversion)
rs_n, sd_n, sd_sigma_n = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                            -0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)
rs_s, sd_s, sd_sigma_s = \
    surfdens_radial_profile(gal, mom0=mom0,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)

rs_feath, sd_feath, sd_sigma_feath = \
    surfdens_radial_profile(gal, mom0=mom0_feath, dr=dr,
                            restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)
rs_feath_n, sd_feath_n, sd_sigma_feath_n = \
    surfdens_radial_profile(gal, mom0=mom0_feath,
                            pa_bounds=Angle([0.5 * np.pi * u.rad,
                                            -0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)
rs_feath_s, sd_feath_s, sd_sigma_feath_s = \
    surfdens_radial_profile(gal, mom0=mom0_feath,
                            pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                             0.5 * np.pi * u.rad]),
                            dr=dr, restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)

# Arecibo
arecibo_file = arecibo_HI_data_path("M33only_jy_stokes_vrad.fits")
arecibo_cube = SpectralCube.read(arecibo_file)

arecibo_mom0 = arecibo_cube.moment0()
rs_arec, sd_arec, sd_sigma_arec = \
    surfdens_radial_profile(gal, cube=arecibo_cube, mom0=arecibo_mom0,
                            dr=dr,
                            restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)

# GBT
gbt_file = gbt_HI_data_path("14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088.fits")
gbt_cube = SpectralCube.read(gbt_file)

gbt_mom0 = gbt_cube.moment0().to(u.K * u.km / u.s)
rs_gbt, sd_gbt, sd_sigma_gbt = \
    surfdens_radial_profile(gal, cube=gbt_cube, mom0=gbt_mom0,
                            dr=dr,
                            restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)
# Archival HI
arch_vla_file = os.path.join(data_path, "VLA/AT0206/old_imaging/m33_hi.masked.fits")
# arch_vla_file = os.path.join(data_path, "VLA/AT0206/imaging/M33_206_b_c_HI.fits")
arch_cube = SpectralCube.read(arch_vla_file)

arch_mom0 = arch_cube.moment0()
rs_arch, sd_arch, sd_sigma_arch = \
    surfdens_radial_profile(gal, cube=arch_cube, mom0=arch_mom0, dr=dr,
                            restfreq=hi_freq,
                            mass_conversion=hi_mass_conversion)

onecolumn_figure(font_scale=1.2)
# Show the total radial profile VLA and Arecibo
p.errorbar(rs.value, sd.value,
           yerr=sd_sigma.value, fmt="-",
           label="VLA", drawstyle='steps-mid')
p.errorbar(rs_feath.value, sd_feath.value,
           yerr=sd_sigma_feath.value, fmt=":",
           label="VLA + GBT", drawstyle='steps-mid')
# p.plot(rs_arec.value, sd_arec.value, "g--", drawstyle='steps-mid',
#        label="Arecibo")
p.plot(rs_gbt.value, sd_gbt.value, "--", drawstyle='steps-mid',
       label="GBT")
# p.errorbar(rs_arec.value, sd_arec.value, yerr=sd_sigma_arec.value,
#            fmt="o--", color="g", label="Arecibo", drawstyle='steps-mid')
p.ylabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid("on")
p.tight_layout()
p.savefig(allfigs_path("HI_properties/M33_surfdens_profile_w_GBT.pdf"))
p.savefig(allfigs_path("HI_properties/M33_surfdens_profile_w_GBT.png"))
p.close()

# W/ archival VLA
ax = p.subplot(111)
p.errorbar(rs.value, sd.value,
           yerr=sd_sigma.value, fmt="-",
           label="VLA", drawstyle='steps-mid')
p.errorbar(rs_feath.value, sd_feath.value,
           yerr=sd_sigma_feath.value, fmt="-.",
           label="VLA + GBT", drawstyle='steps-mid')
p.plot(rs_gbt.value, sd_gbt.value, "--", drawstyle='steps-mid',
       label="GBT")
# p.plot(rs_arec.value, sd_arec.value, "g--", drawstyle='steps-mid',
#        label="Arecibo")
p.errorbar(rs_arch.value, sd_arch.value, yerr=sd_sigma_arec.value,
           fmt=":", drawstyle='steps-mid',
           label="Archival VLA")

ax.add_patch(patches.Rectangle((1.9, 0.3), 2.9, 1.4, facecolor='w',
                               edgecolor='k'))
# Add lines according to the beam widths
conv = 4.e-3 * u.kpc / u.arcsec
p.plot([2, 2 + (gbt_cube.beam.major.to(u.arcsec) * conv).value],
       [1.5, 1.5], color=sb.color_palette()[2])
# p.plot([2, 2 + (arec_cube.beam.major.to(u.arcsec) * conv).value],
#        [1.5, 1.5], color=sb.color_palette()[2])
p.plot([2, 2 + (mom0.beam.major.to(u.arcsec) * conv).value],
       [1.0, 1.0], color=sb.color_palette()[0])
p.plot([2, 2 + (arch_cube.beam.major.to(u.arcsec) * conv).value],
       [0.5, 0.5], color=sb.color_palette()[3])
p.ylabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.ylim([0, 10])
p.xlim([0, 13])
p.legend(loc='upper right', frameon=True)
p.grid("on")
p.tight_layout()

p.savefig(allfigs_path("HI_properties/M33_Sigma_profile_w_GBT_archival.pdf"))
p.savefig(allfigs_path("HI_properties/M33_Sigma_profile_w_GBT_archival.png"))
p.close()

# Show the north vs south profiles
cpal = sb.color_palette()

p.plot(rs.value, sd.value, color='k',
       drawstyle='steps-mid', label="Total")
p.plot(rs_n.value, sd_n.value, "-.", label="North",
       drawstyle='steps-mid', color=cpal[2])
p.plot(rs_s.value, sd_s.value, "--", label="South",
       drawstyle='steps-mid', color=cpal[0])
# p.errorbar(rs_n.value, sd_n.value, yerr=sd_sigma_n.value, fmt="D-",
#            color="b", label="North")
# p.errorbar(rs_s.value, sd_s.value, yerr=sd_sigma_s.value, fmt="o-",
#            color="g", label="South")
p.ylabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid("on")
p.ylim([0, 8])
p.xlim([0, 10])
p.tight_layout()

p.savefig(allfigs_path("HI_properties/M33_Sigma_profile_N_S.pdf"))
p.savefig(allfigs_path("HI_properties/M33_Sigma_profile_N_S.png"))

p.close()

p.plot(rs_feath.value, sd_feath.value, color='k',
       drawstyle='steps-mid', label="Total")
p.plot(rs_feath_n.value, sd_feath_n.value, "-.", label="North",
       drawstyle='steps-mid', color=cpal[2])
p.plot(rs_feath_s.value, sd_feath_s.value, "--", label="South",
       drawstyle='steps-mid', color=cpal[0])
p.ylabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid("on")
p.ylim([0, 10.5])
p.xlim([0, 10])
p.fill_between([0, 0.5], 0, 11, color='gray', alpha=0.5)
p.tight_layout()
p.savefig(allfigs_path("HI_properties/M33_feathered_Sigma_profile_N_S.pdf"))
p.savefig(allfigs_path("HI_properties/M33_feathered_Sigma_profile_N_S.png"))
p.close()

# Compare to the surface density profile in Corbelli
corbelli_filename = \
    c_hi_analysispath("rotation_curves/corbelli_rotation_curve.csv")
corbelli = Table.read(corbelli_filename)

p.plot(rs.value, sd.value,
       linestyle="-",
       label="This work - VLA", drawstyle='steps-mid')
p.plot(rs_feath.value, sd_feath.value,
       linestyle=":",
       label="This work - VLA + GBT", drawstyle='steps-mid')
p.plot(corbelli["R"][corbelli["R"] <= 10.0],
       corbelli["SigmaHI"][corbelli["R"] <= 10.0], "--",
       drawstyle='steps-mid',
       label="Corbelli et al. (2014)")
p.ylabel(r"$\Sigma$ (M$_{\odot}$ pc$^{-2}$)")
p.xlabel(r"Radius (kpc)")
p.legend(loc='best', frameon=True)
p.grid()
p.tight_layout()
p.savefig(allfigs_path("HI_properties/M33_Sigma_profile_w_Corbelli.pdf"))
p.savefig(allfigs_path("HI_properties/M33_Sigma_profile_w_Corbelli.png"))
p.close()

default_figure()
