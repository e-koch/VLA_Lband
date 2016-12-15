
'''
Compare rotation curves found from the centroid, peak velocity, and GH peak
surfaces.
'''

import matplotlib.pyplot as plt
import astropy.wcs as wcs
from astropy.io import fits
import astropy.units as u
import pandas as pd
import numpy as np

from paths import fourteenB_HI_data_path, paper1_figures_path
from constants import pb_lim
from galaxy_params import gal


# Requires that diskfit has already been run on these 3 surfaces!
# By default, the analysis pipeline should run each.
cent_name = "diskfit_noasymm_noradial_nowarp_output/rad.out.csv"
peak_name = "diskfit_peakvels_noasymm_noradial_nowarp_output/rad.out.csv"
gh_name = "diskfit_ghfit_noasymm_noradial_nowarp_output/rad.out.csv"

cent_tab = pd.read_csv(fourteenB_HI_data_path(cent_name))
peak_tab = pd.read_csv(fourteenB_HI_data_path(peak_name))
peak_tab["Vt"] = np.abs(peak_tab["Vt"])
gh_tab = pd.read_csv(fourteenB_HI_data_path(gh_name))

# Read in WCS info to convert to distances
mom1_name = fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_{}.ellip_mask.mom1.fits".format(pb_lim))
mom1 = fits.open(mom1_name)
mom1_wcs = wcs.WCS(mom1[0].header)

scale = wcs.utils.proj_plane_pixel_scales(mom1_wcs)[0]
# Distance scaling (1" ~ 4 pc). Conversion is deg to kpc
dist_scale = (np.pi / 180.) * gal.distance.to(u.kpc).value


# Plot the rotation curves
plt.errorbar(cent_tab["r"] * scale * dist_scale, cent_tab["Vt"],
             yerr=cent_tab["eVt"], markersize=7,
             label='Centroid', fmt='D-', alpha=0.8)

plt.errorbar(peak_tab["r"] * scale * dist_scale, peak_tab["Vt"],
             yerr=peak_tab["eVt"], markersize=7,
             label='Peak Velocity', fmt='o-', alpha=0.8)

plt.errorbar(gh_tab["r"] * scale * dist_scale, gh_tab["Vt"],
             yerr=gh_tab["eVt"], markersize=7,
             label='Gauss-Hermite Fit', fmt='+-', alpha=0.8)

plt.legend(loc='lower right')
plt.grid()

plt.ylabel(r"V$_{\rm circ}$ (km/s)")
plt.xlabel("Radius (kpc)")

plt.savefig(paper1_figures_path("HI_rotation_curve_comparisons.png"))
plt.savefig(paper1_figures_path("HI_rotation_curve_comparisons.pdf"))

plt.close()
