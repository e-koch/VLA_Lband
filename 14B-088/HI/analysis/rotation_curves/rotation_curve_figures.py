
'''
Make figures from the outputs from fit_rotation_curves.py and
fit_rotation_curves_feathered.py
'''

import matplotlib.pyplot as plt
from galaxies import Galaxy
from astropy.table import Table
from astropy import wcs
from astropy.io import fits
import astropy.units as u

from cube_analysis.rotation_curves import generate_vrot_model
from cube_analysis.rotation_curves.curve_fitting import vcirc

from paths import (fourteenB_HI_data_path, fourteenB_HI_data_wGBT_path,
                   fourteenB_HI_file_dict, fourteenB_wGBT_HI_file_dict,
                   paper1_figures_path, allfigs_path,
                   c_hi_analysispath)
from constants import ang_to_phys

from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_figure, twocolumn_twopanel_figure)


gal = Galaxy("M33")

# Need the grid size
mom1_data = fits.open(fourteenB_HI_file_dict['Moment1'])[0]
scale = wcs.utils.proj_plane_pixel_scales(wcs.WCS(mom1_data.header))[0] * u.deg

# Load in the tables
mom1_vlaonly = fourteenB_HI_data_path("diskfit_noasymm_noradial_nowarp_output/rad.out.csv")
peakvel_vlaonly = fourteenB_HI_data_path("diskfit_peakvels_noasymm_noradial_nowarp_output/rad.out.csv")

mom1_feath = fourteenB_HI_data_wGBT_path("diskfit_noasymm_noradial_nowarp_output/rad.out.csv")
peakvel_feath = fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_noradial_nowarp_output/rad.out.csv")

print(argh)
# onecolumn_figure(font_scale=1.0)
twocolumn_figure()

# Four panel figure of the different models
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)

labels = ["Centroid VLA-only", "Peak Vel. VLA-only",
          "Centroid VLA + GBT", "Peak Vel. VLA + GBT"]
tables = [mom1_vlaonly, peakvel_vlaonly, mom1_feath, peakvel_feath]

i = 0
for ax, fname, label in zip(axes.ravel(), tables, labels):

    data = Table.read(fname)

    pars, pcov = generate_vrot_model(data)

    phys_radius = ang_to_phys(data["r"] * scale).to(u.kpc).value
    plot_pars = pars.copy()
    plot_pars[2] *= ang_to_phys(scale).to(u.kpc).value
    ax.errorbar(phys_radius, data["Vt"], yerr=data['eVt'],
                fmt='-', label="DiskFit",
                drawstyle='steps-mid')
    ax.plot(phys_radius, vcirc(phys_radius, *plot_pars), '-', linewidth=2,
            label='Model')
    ax.text(4, 40, label, bbox={"boxstyle": "square", "facecolor": "w"})
    if i == 0 or i == 2:
        ax.set_ylabel(r"Circular Velocity (km / s)")
    if i == 2 or i == 3:
        ax.set_xlabel(r"Radius (kpc)")
    ax.grid()

    if i == 0:
        ax.legend(loc='upper left', frameon=True)

    i += 1

plt.tight_layout()

plt.savefig(allfigs_path("M33_vrot_models_wfit.pdf"))
plt.savefig(allfigs_path("M33_vrot_models_wfit.png"))
plt.close()

# Plot all four rotation curves in one panel, and then they're fits in the
# other
twocolumn_twopanel_figure()
fig, ax = plt.subplots(1, 2, sharey=True, sharex=True)
symbols = ["-", "--"]
for fname, label, symbol in zip(tables, labels, symbols * 2):
    data = Table.read(fname)

    pars, pcov = generate_vrot_model(data)

    phys_radius = ang_to_phys(data["r"] * scale).to(u.kpc).value
    plot_pars = pars.copy()
    plot_pars[2] *= ang_to_phys(scale).to(u.kpc).value
    ax[0].errorbar(phys_radius, data["Vt"], yerr=data['eVt'],
                   fmt=symbol, label=label,
                   drawstyle='steps-mid', alpha=1.0)

    ax[1].plot(phys_radius, vcirc(phys_radius, *plot_pars), symbol,
               linewidth=2, alpha=1.0)
ax[0].grid()
ax[1].grid()
ax[0].set_ylabel(r"Circular Velocity (km / s)")
ax[0].set_xlabel(r"Radius (kpc)")
ax[1].set_xlabel(r"Radius (kpc)")
ax[0].legend(loc='lower right', frameon=True)
plt.tight_layout()

plt.savefig(allfigs_path("M33_vrot_models_wfit_comparisons.pdf"))
plt.savefig(allfigs_path("M33_vrot_models_wfit_comparisons.png"))
plt.close()

# Use the peak velocity model from the feathered cube. It's essentially
# identical to the VLA-only peak velocity curve.

# load in the Corbelli and Kam curves for comparison
onecolumn_figure()
corbelli = Table.read(c_hi_analysispath("rotation_curves/corbelli_rotation_curve.csv"))
kam = Table.read(c_hi_analysispath("rotation_curves/kam_rotation_curve.csv"))

data = Table.read(peakvel_feath)

pars, pcov = generate_vrot_model(data)

phys_radius = ang_to_phys(data["r"] * scale).to(u.kpc).value
plot_pars = pars.copy()
plot_pars[2] *= ang_to_phys(scale).to(u.kpc).value
plt.errorbar(phys_radius, data["Vt"], yerr=data['eVt'],
             fmt="-", label="This Work",
             drawstyle='steps-mid')
plt.errorbar(corbelli["R"][corbelli["R"] <= 8.0],
             corbelli["Vr"][corbelli["R"] <= 8.0],
             yerr=corbelli["dVr"][corbelli["R"] <= 8.0],
             fmt='--', label="Corbelli et al. (2014)",
             drawstyle='steps-mid')
plt.errorbar(kam["R"][kam["R"] <= 8.0],
             kam["Vr"][kam["R"] <= 8.0],
             yerr=kam["dVr"][kam["R"] <= 8.0],
             fmt=':', label="Kam et al. (2017)",
             drawstyle='steps-mid')
plt.ylabel(r"Circular Velocity (km / s)")
plt.xlabel(r"Radius (kpc)")
plt.legend(loc='lower right', frameon=True)
plt.grid()

plt.tight_layout()

plt.savefig(paper1_figures_path("M33_vrot_feathered_peakvels_wCorbelli_Kam.pdf"))
plt.savefig(paper1_figures_path("M33_vrot_feathered_peakvels_wCorbelli_Kam.png"))

plt.close()

default_figure()
