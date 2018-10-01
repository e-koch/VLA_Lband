
'''
Create stacked profiles based on the clean sample of per-los
Gaussian fitting.

- Stacked profiles of HI and CO based on the fitted components, centered
  on the HI
- Same, but centered on the CO peak
- Stacked profiles of HI and CO using the actual spectra from the fitting,
  centered on the HI
'''

from spectral_cube import SpectralCube
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from pandas import DataFrame
import astropy.units as u
from astropy.utils.console import ProgressBar
import os
from os.path import join as osjoin
import seaborn as sb
from scipy.special import erf

from cube_analysis.spectral_stacking_models import find_hwhm, fit_gaussian, fit_hwhm
from cube_analysis.spectra_shifter import cube_shifter
from cube_analysis.spectral_fitting import gauss_model_discrete

from paths import (fourteenB_wGBT_HI_file_dict, iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path, allfigs_path)
from plotting_styles import (default_figure, onecolumn_figure,
                             onecolumn_twopanel_figure)
from galaxy_params import gal_feath as gal
from constants import (co21_mass_conversion, hi_mass_conversion, hi_freq,
                       beam_eff_30m_druard)


default_figure()

tab = Table.read(fourteenB_HI_data_wGBT_path("tables/hi_co_gaussfit_column_densities_perpix.fits"))

good_pts = np.logical_and(~tab['multicomp_flag_HI'],
                          ~tab['multicomp_flag_CO'])
good_pts = np.logical_and(good_pts,
                          tab["sigma_HI"] > 3800)
# Minimum CO line width of one channel.
good_pts = np.logical_and(good_pts,
                          tab["sigma_CO"] >= 2600)


cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])
co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.fits"))


co_mask = fits.open(iram_co21_14B088_data_path(
    "m33.co21_iram.14B-088_HI_source_mask.fits"))[0].data

# Require there be 3 pixels in the mask
co_mask_spatial = co_mask.sum(0) > 3
yposns, xposns = np.where(co_mask_spatial)

hi_specaxis = cube.spectral_axis
co_specaxis = co_cube.spectral_axis

hi_chanwidth = np.abs(np.diff(hi_specaxis[:2])[0].value)
co_chanwidth = np.abs(np.diff(co_specaxis[:2])[0].value)

hi_mod_vals = {}
co_mod_vals = {}

# Create stacked profiles from the fitted parameters

hi_stack_modelspec_wHI = np.zeros_like(cube[:, 0, 0].value)
co_stack_modelspec_wHI = np.zeros_like(co_cube[:, 0, 0].value)

# Loop through model parameters

# The HI centre is fixed in this case
hi_centre = cube.spectral_axis[hi_stack_modelspec_wHI.size // 2].value

for i in ProgressBar(np.where(good_pts)[0]):

    # The CO centre is the fitted centre relative to the HI
    co_centre = (tab['mean_CO'][i] - tab['mean_HI'][i]) + hi_centre

    hi_stack_modelspec_wHI += \
        gauss_model_discrete(cube.spectral_axis.value,
                             tab['amp_HI'][i], hi_centre, tab['sigma_HI'][i])

    co_stack_modelspec_wHI += \
        gauss_model_discrete(co_cube.spectral_axis.value,
                             tab['amp_CO'][i], co_centre, tab['sigma_CO'][i])


# Fit the HWHM model to these

hi_model_hwhm_wHI = fit_hwhm(cube.spectral_axis.value, hi_stack_modelspec_wHI,
                             niters=1000)

hi_mod_vals["model_HImodelvels"] = hi_model_hwhm_wHI[0]
hi_mod_vals["model_HImodelvels_lower"] = hi_model_hwhm_wHI[1][0]
hi_mod_vals["model_HImodelvels_upper"] = hi_model_hwhm_wHI[1][1]

co_model_hwhm_wHI = fit_hwhm(co_cube.spectral_axis.value,
                             co_stack_modelspec_wHI,
                             niters=1000)

co_mod_vals["model_HImodelvels"] = co_model_hwhm_wHI[0]
co_mod_vals["model_HImodelvels_lower"] = co_model_hwhm_wHI[1][0]
co_mod_vals["model_HImodelvels_upper"] = co_model_hwhm_wHI[1][1]

# Do the same thing, now treating CO as the line centre

hi_stack_modelspec_wCO = np.zeros_like(cube[:, 0, 0].value)
co_stack_modelspec_wCO = np.zeros_like(co_cube[:, 0, 0].value)

# Loop through model parameters

# The HI centre is fixed in this case
co_centre = co_cube.spectral_axis[co_stack_modelspec_wCO.size // 2].value

for i in ProgressBar(np.where(good_pts)[0]):

    # The CO centre is the fitted centre relative to the HI
    hi_centre = (tab['mean_HI'][i] - tab['mean_CO'][i]) + co_centre

    hi_stack_modelspec_wCO += \
        gauss_model_discrete(cube.spectral_axis.value,
                             tab['amp_HI'][i], hi_centre, tab['sigma_HI'][i])

    co_stack_modelspec_wCO += \
        gauss_model_discrete(co_cube.spectral_axis.value,
                             tab['amp_CO'][i], co_centre, tab['sigma_CO'][i])


# Fit the HWHM model to these

hi_model_hwhm_wCO = fit_hwhm(cube.spectral_axis.value, hi_stack_modelspec_wCO,
                             niters=1000)

hi_mod_vals["model_COmodelvels"] = hi_model_hwhm_wCO[0]
hi_mod_vals["model_COmodelvels_lower"] = hi_model_hwhm_wCO[1][0]
hi_mod_vals["model_COmodelvels_upper"] = hi_model_hwhm_wCO[1][1]


co_model_hwhm_wCO = fit_hwhm(co_cube.spectral_axis.value,
                             co_stack_modelspec_wCO,
                             niters=1000)

co_mod_vals["model_COmodelvels"] = co_model_hwhm_wCO[0]
co_mod_vals["model_COmodelvels_lower"] = co_model_hwhm_wCO[1][0]
co_mod_vals["model_COmodelvels_upper"] = co_model_hwhm_wCO[1][1]


# Now perform stacking using the actual spectra, and not the models

print("Stacking full spectra in good sample")

yposns = tab['ypts'][good_pts]
xposns = tab['xpts'][good_pts]

peakvels_fit_HI = np.zeros_like(cube[0].value)

peakvels_fit_HI[yposns, xposns] = tab['mean_HI'][good_pts]
peakvels_fit_HI[peakvels_fit_HI == 0.] = np.NaN


shifted_spectra_HI = cube_shifter(cube, peakvels_fit_HI * u.m / u.s,
                                  xy_posns=(yposns, xposns),
                                  return_spectra=True,
                                  verbose=True)

shifted_spectra_CO = cube_shifter(co_cube, peakvels_fit_HI * u.m / u.s,
                                  xy_posns=(yposns, xposns),
                                  return_spectra=True,
                                  verbose=True)

hi_stack_spec_wHI = sum(shifted_spectra_HI[0])
co_stack_spec_wHI = sum(shifted_spectra_CO[0])

# Fit the HWHM model to these

hi_model_hwhm_spec_wHI = fit_hwhm(hi_stack_spec_wHI.spectral_axis.value,
                                  hi_stack_spec_wHI.value,
                                  niters=1000)

hi_mod_vals["spec_HImodelvels"] = hi_model_hwhm_spec_wHI[0]
hi_mod_vals["spec_HImodelvels_lower"] = hi_model_hwhm_spec_wHI[1][0]
hi_mod_vals["spec_HImodelvels_upper"] = hi_model_hwhm_spec_wHI[1][1]


co_model_hwhm_spec_wHI = fit_hwhm(co_stack_spec_wHI.spectral_axis.value,
                                  co_stack_spec_wHI.value,
                                  niters=1000)

co_mod_vals["spec_HImodelvels"] = co_model_hwhm_spec_wHI[0]
co_mod_vals["spec_HImodelvels_lower"] = co_model_hwhm_spec_wHI[1][0]
co_mod_vals["spec_HImodelvels_upper"] = co_model_hwhm_spec_wHI[1][1]


# And stacking actual spectra with the measured (not fit) peak velocity

print("Stacking full spectra in good sample w/ original HI peak velocities")

hi_peakvels = fits.open(fourteenB_wGBT_HI_file_dict['PeakVels'])[0].data / 1000.

yposns = tab['ypts'][good_pts]
xposns = tab['xpts'][good_pts]

peakvels_LOS_HI = np.zeros_like(cube[0].value)

peakvels_LOS_HI[yposns, xposns] = hi_peakvels[yposns, xposns]
peakvels_LOS_HI[peakvels_fit_HI == 0.] = np.NaN


shifted_spectra_orgvelcent_HI = cube_shifter(cube, peakvels_LOS_HI * u.km / u.s,
                                             xy_posns=(yposns, xposns),
                                             return_spectra=True,
                                             verbose=True)

shifted_spectra_orgvelcent_CO = cube_shifter(co_cube, peakvels_LOS_HI * u.km / u.s,
                                             xy_posns=(yposns, xposns),
                                             return_spectra=True,
                                             verbose=True)

hi_stack_all_spec_wHI = sum(shifted_spectra_orgvelcent_HI[0])
co_stack_all_spec_wHI = sum(shifted_spectra_orgvelcent_CO[0])

# Fit the HWHM model to these

hi_model_hwhm_orgvelcent_wHI = fit_hwhm(hi_stack_all_spec_wHI.spectral_axis.value,
                                        hi_stack_all_spec_wHI.value,
                                        niters=1000)

hi_mod_vals["spec_HIpeakvels"] = hi_model_hwhm_orgvelcent_wHI[0]
hi_mod_vals["spec_HIpeakvels_lower"] = hi_model_hwhm_orgvelcent_wHI[1][0]
hi_mod_vals["spec_HIpeakvels_upper"] = hi_model_hwhm_orgvelcent_wHI[1][1]

co_model_hwhm_orgvelcent_wHI = \
    fit_hwhm(co_stack_all_spec_wHI.spectral_axis.value,
             co_stack_all_spec_wHI.value,
             niters=1000)

co_mod_vals["spec_HIpeakvels"] = co_model_hwhm_orgvelcent_wHI[0]
co_mod_vals["spec_HIpeakvels_lower"] = co_model_hwhm_orgvelcent_wHI[1][0]
co_mod_vals["spec_HIpeakvels_upper"] = co_model_hwhm_orgvelcent_wHI[1][1]

# Compare stacking with all LOS in the set
# Poor results due to bad fits. Skip this comparison

# print("Stacking full spectra in total sample")

# yposns = tab['ypts']
# xposns = tab['xpts']

# peakvels_fit_HI = np.zeros_like(cube[0].value)

# peakvels_fit_HI[yposns, xposns] = tab['mean_HI']
# peakvels_fit_HI[peakvels_fit_HI == 0.] = np.NaN


# shifted_spectra_all_HI = cube_shifter(cube, peakvels_fit_HI * u.m / u.s,
#                                       xy_posns=(yposns, xposns),
#                                       return_spectra=True,
#                                       verbose=True)

# shifted_spectra_all_CO = cube_shifter(co_cube, peakvels_fit_HI * u.m / u.s,
#                                       xy_posns=(yposns, xposns),
#                                       return_spectra=True,
#                                       verbose=True)

# hi_stack_all_spec_wHI = sum(shifted_spectra_all_HI[0])
# co_stack_all_spec_wHI = sum(shifted_spectra_all_CO[0])

# # Fit the HWHM model to these

# hi_model_hwhm_all_spec_wHI = fit_hwhm(hi_stack_all_spec_wHI.spectral_axis.value,
#                                       hi_stack_all_spec_wHI.value,
#                                       niters=1000)

# hi_mod_vals["spec_HImodelvels_fullsamp"] = hi_model_hwhm_all_spec_wHI[0]
# hi_mod_vals["spec_HImodelvels_fullsamp_lower"] = hi_model_hwhm_all_spec_wHI[1][0]
# hi_mod_vals["spec_HImodelvels_fullsamp_upper"] = hi_model_hwhm_all_spec_wHI[1][1]


# co_model_hwhm_all_spec_wHI = fit_hwhm(co_stack_all_spec_wHI.spectral_axis.value,
#                                       co_stack_all_spec_wHI.value,
#                                       niters=1000)

# co_mod_vals["spec_HImodelvels_fullsamp"] = co_model_hwhm_all_spec_wHI[0]
# co_mod_vals["spec_HImodelvels_fullsamp_lower"] = co_model_hwhm_all_spec_wHI[1][0]
# co_mod_vals["spec_HImodelvels_fullsamp_upper"] = co_model_hwhm_all_spec_wHI[1][1]


# And stacking actual spectra with the measured (not fit) peak velocity

print("Stacking full spectra in total sample w/ original HI peak velocities")

hi_peakvels = fits.open(fourteenB_wGBT_HI_file_dict['PeakVels'])[0].data / 1000.

yposns = tab['ypts']
xposns = tab['xpts']

peakvels_LOS_HI = np.zeros_like(cube[0].value)

peakvels_LOS_HI[yposns, xposns] = hi_peakvels[yposns, xposns]
peakvels_LOS_HI[peakvels_fit_HI == 0.] = np.NaN


shifted_spectra_all_orgvelcent_HI = cube_shifter(cube,
                                                 peakvels_LOS_HI * u.km / u.s,
                                                 xy_posns=(yposns, xposns),
                                                 return_spectra=True,
                                                 verbose=True)

shifted_spectra_all_orgvelcent_CO = cube_shifter(co_cube,
                                                 peakvels_LOS_HI * u.km / u.s,
                                                 xy_posns=(yposns, xposns),
                                                 return_spectra=True,
                                                 verbose=True)

hi_stack_all_spec_orgvelcent_wHI = sum(shifted_spectra_all_orgvelcent_HI[0])
co_stack_all_spec_orgvelcent_wHI = sum(shifted_spectra_all_orgvelcent_CO[0])

# Fit the HWHM model to these

hi_model_hwhm_all_spec_orgvelcent_wHI = \
    fit_hwhm(hi_stack_all_spec_orgvelcent_wHI.spectral_axis.value,
             hi_stack_all_spec_orgvelcent_wHI.value,
             niters=1000)

hi_mod_vals["spec_HIpeakvels_fullsamp"] = hi_model_hwhm_all_spec_orgvelcent_wHI[0]
hi_mod_vals["spec_HIpeakvels_fullsamp_lower"] = hi_model_hwhm_all_spec_orgvelcent_wHI[1][0]
hi_mod_vals["spec_HIpeakvels_fullsamp_upper"] = hi_model_hwhm_all_spec_orgvelcent_wHI[1][1]

co_model_hwhm_all_spec_orgvelcent_wHI = \
    fit_hwhm(co_stack_all_spec_orgvelcent_wHI.spectral_axis.value,
             co_stack_all_spec_orgvelcent_wHI.value,
             niters=1000)

co_mod_vals["spec_HIpeakvels_fullsamp"] = co_model_hwhm_all_spec_orgvelcent_wHI[0]
co_mod_vals["spec_HIpeakvels_fullsamp_lower"] = co_model_hwhm_all_spec_orgvelcent_wHI[1][0]
co_mod_vals["spec_HIpeakvels_fullsamp_upper"] = co_model_hwhm_all_spec_orgvelcent_wHI[1][1]

# Save these stacking test parameters in a nice table

parnames_hwhm = co_model_hwhm_all_spec_orgvelcent_wHI[2]

df_hi = DataFrame(hi_mod_vals, index=parnames_hwhm)
# df_hi.to_latex(alltables_path("hi_totalprof_hwhm_lossample_stacking.tex"))
df_hi.to_csv(fourteenB_HI_data_wGBT_path("tables/hi_totalprof_hwhm_lossample_stacking.csv",
                                         no_check=True))

df_co = DataFrame(co_mod_vals, index=parnames_hwhm)
# df_co.to_latex(alltables_path("co_totalprof_hwhm_lossample_stacking.tex"))
df_co.to_csv(iram_co21_14B088_data_path("tables/co_totalprof_hwhm_lossample_stacking.csv",
                                        no_check=True))

# Compare the different stacking approaches

figure_folder = allfigs_path("stacked_profiles")
if not os.path.exists(figure_folder):
    os.mkdir(figure_folder)

# onecolumn_figure()
default_figure()

val_names = ['model_COmodelvels', 'model_HImodelvels', 'spec_HImodelvels',
             # 'spec_HImodelvels_fullsamp',
             'spec_HIpeakvels',
             'spec_HIpeakvels_fullsamp']

low_names = ["{}_lower" for val in val_names]
up_names = ["{}_upper" for val in val_names]

fig, axs = plt.subplots(2, 2, sharex=True)

plot_cols = ['sigma', 'f_wings', 'asymm', 'kappa']
label_cols = ['sigma', r'$f_\mathrm{wings}$', 'asymm', 'kappa']

for ax, col, lab in zip(axs.ravel(), plot_cols, label_cols):
    ax.errorbar(range(5), df_hi.T[col][val_names],
                yerr=[df_hi.T[col][low_names],
                      df_hi.T[col][up_names]], label='HI')
    ax.errorbar(range(5), df_co.T[col][val_names],
                yerr=[df_co.T[col][low_names],
                      df_co.T[col][up_names]], label='CO')
    ax.set_ylabel(lab)
    ax.grid()

axs.ravel()[0].legend(frameon=True)

axs.ravel()[-1].set_xticklabels([val.replace("_", " ") for val in val_names],
                                rotation=45)
axs.ravel()[-2].set_xticklabels([val.replace("_", " ") for val in val_names],
                                rotation=45)

plt.tight_layout()

fig.savefig(osjoin(figure_folder, "HI_CO_stacked_profile_LOSmodel_comparisons.pdf"))
fig.savefig(osjoin(figure_folder, "HI_CO_stacked_profile_LOSmodel_comparisons.png"))
plt.close()
