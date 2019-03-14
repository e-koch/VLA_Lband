
'''
Compare pressure models on GMC-scales in M33.
'''


from astropy.io import fits
from spectral_cube import Projection
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import seaborn as sb
from corner import hist2d
from os.path import join as osjoin
import os
import matplotlib.lines as mlines
from astropy.table import Table
from scipy import odr
from astropy.visualization import hist as astro_hist
from astropy.stats import histogram
from astropy.table import Table

from paths import (data_path, fourteenB_wGBT_HI_file_dict,
                   fourteenB_HI_data_wGBT_path,
                   iram_co21_14B088_data_path, allfigs_path)
# from paths import data_path
from constants import (co21_mass_conversion, hi_mass_conversion, hi_freq,
                       beam_eff_30m_druard)
from plotting_styles import (default_figure, onecolumn_figure,
                             twocolumn_figure,
                             onecolumn_Npanel_figure,
                             twocolumn_twopanel_figure,
                             fullpage_figure)
from galaxy_params import gal_feath as gal


from cube_analysis.h2_models import (Pext_star, Pext_sum,
                                     krumholz2013_sigmaHI_H2,
                                     sternberg2014_model_const_alphG,
                                     sternberg2014_model_const_alphG_avg,
                                     sternberg2014_sigmaHI_H2,
                                     )
from cube_analysis.eiv_emcee import bayes_linear


fig_path = osjoin(allfigs_path(""), "co_vs_hi/h2_hi_stellar_pressure")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)

# Load radial table from radial_surfdens_avg_table.py
tab = Table.read(fourteenB_HI_data_wGBT_path("tables/sd_radial_500pc_w_stellar_dust.fits"))

dust_col_hdu = fits.open(osjoin(data_path, "m33_dust.surface.density_FB.beta=1.8_gauss41.0_regrid_bksub.fits"))[0]
dust_temp_hdu = fits.open(osjoin(data_path, "m33_dust.temperature_FB.beta=1.8_gauss41.0_regrid_bksub.fits"))[0]

hi_mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0']))
co_mom0 = Projection.from_hdu(fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits")))

# In K m s-1
hi_mom0_sigma = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.moment0_Kms_err.fits")))
hi_mom0_sigma = hi_mom0_sigma.to(u.K * u.km / u.s) * np.cos(gal.inclination)


# Convert to K km/s in HI
hi_mom0_Kkms = (hi_mom0 * hi_mom0.beam.jtok(hi_freq) * u.beam / u.Jy).to(u.K * u.km / u.s) * np.cos(gal.inclination)

hi_surfdens = hi_mass_conversion * hi_mom0_Kkms
hi_surfdens_sigma = hi_mass_conversion * hi_mom0_sigma

co_mom0 = co_mom0.to(u.K * u.km / u.s)
co_surfdens = co_mom0 * co21_mass_conversion * np.cos(gal.inclination) / beam_eff_30m_druard

co_mask = fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_source_mask.fits"))[0]

# good_pts = np.where(np.isfinite(mom0_reproj))
good_co_mask = co_mask.data.sum(0) >= 2

co_surfdens = co_surfdens.with_mask(good_co_mask)

# Radius
radius = gal.radius(header=hi_mom0.header).to(u.kpc)

dust_col_hdr = dust_col_hdu.header.copy()
dust_col_hdr['BUNIT'] = "solMass / pc2"
dust_surfdens = Projection.from_hdu(fits.PrimaryHDU(dust_col_hdu.data[0],
                                                    dust_col_hdr))
dust_surfdens_rep = dust_surfdens.reproject(hi_mom0.header)

dust_surfdens_sigma = Projection.from_hdu(fits.PrimaryHDU(dust_col_hdu.data[2],
                                                          dust_col_hdr))
dust_surfdens_sigma_rep = dust_surfdens_sigma.reproject(hi_mom0.header)

# Apply inclination correction
dust_surfdens_rep = dust_surfdens_rep * np.cos(gal.inclination)
dust_surfdens_sigma_rep = dust_surfdens_sigma_rep * np.cos(gal.inclination)


dust_temp_hdr = dust_temp_hdu.header.copy()
dust_temp_hdr['BUNIT'] = "K"
dust_temp = Projection.from_hdu(fits.PrimaryHDU(dust_temp_hdu.data[0],
                                                dust_temp_hdr))
dust_temp_rep = dust_temp.reproject(hi_mom0.header)

stellar_surfdens = Projection.from_hdu(fits.open(osjoin(data_path, "m33.stellarmass.fits")))
stellar_surfdens = stellar_surfdens * u.solMass / u.kpc**2
stellar_surfdens_rep = stellar_surfdens.reproject(hi_mom0.header).to(u.solMass / u.pc**2)


# Pressure comparisons

co_mask = np.isfinite(co_surfdens) & (co_surfdens.value > 0.5)

hi_mask = np.isfinite(hi_surfdens) & (hi_surfdens.value > 0)


pext_on_gmc = Pext_sum(co_surfdens,
                       hi_surfdens,
                       stellar_surfdens_rep,
                       H_mol=100 * u.pc,
                       H_star=300 * u.pc,
                       sigma_HI=7 * u.km / u.s)

onecolumn_figure()

_ = astro_hist(np.log10(pext_on_gmc[0][co_mask].value), bins='knuth', alpha=0.5,
               label='Total')

_ = astro_hist(np.log10(pext_on_gmc[1][co_mask].value), bins='knuth', alpha=0.5,
               label='GMC')

_ = astro_hist(np.log10(pext_on_gmc[2][co_mask].value), bins='knuth', alpha=0.5,
               label='Stellar on GMC')

_ = astro_hist(np.log10(pext_on_gmc[3][co_mask].value), bins='knuth', alpha=0.5,
               label='HI on GMC')

_ = astro_hist(np.log10(pext_on_gmc[4][co_mask].value), bins='knuth', alpha=0.5,
               label='Stellar on HI')

plt.legend(loc='upper left')
plt.xlabel(r"log Pressure (K cm$^{-3}$)")

plt.savefig(osjoin(fig_path, "pext_total_comparison_hist.pdf"))
plt.savefig(osjoin(fig_path, "pext_total_comparison_hist.png"))
plt.close()

# Do these relations change if we select on Sigma_H2 > Sigma_HI?

plt.figure()

mol_dom_mask = (co_surfdens > hi_surfdens) & co_mask

_ = astro_hist(np.log10(pext_on_gmc[0][mol_dom_mask].value), bins='knuth', alpha=0.5,
               label='Total')

_ = astro_hist(np.log10(pext_on_gmc[1][mol_dom_mask].value), bins='knuth', alpha=0.5,
               label='GMC')

_ = astro_hist(np.log10(pext_on_gmc[2][mol_dom_mask].value), bins='knuth', alpha=0.5,
               label='Stellar on GMC')

_ = astro_hist(np.log10(pext_on_gmc[3][mol_dom_mask].value), bins='knuth', alpha=0.5,
               label='HI on GMC')

_ = astro_hist(np.log10(pext_on_gmc[4][mol_dom_mask].value), bins='knuth', alpha=0.5,
               label='Stellar on HI')

plt.legend(loc='upper left')
plt.xlabel(r"log Pressure (K cm$^{-3}$)")

plt.savefig(osjoin(fig_path, "pext_total_comparison_moldom_hist.pdf"))
plt.savefig(osjoin(fig_path, "pext_total_comparison_moldom_hist.png"))
plt.close()


# The pressure from HI is consistently high here.

# What happens to the HI pressure if the resolution is much lower, compared to
# the same CO resolution?

# Convert to K km/s in HI
hi_mom0_2beam = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("smooth_2beam/M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.mom0.fits")))
hi_mom0_2beam[np.isnan(hi_mom0)] = np.NaN

hi_mom0_5beam = Projection.from_hdu(fits.open(fourteenB_HI_data_wGBT_path("smooth_5beam/M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.mom0.fits")))
hi_mom0_5beam[np.isnan(hi_mom0)] = np.NaN


hi_mom0_2beam_Kkms = (hi_mom0_2beam * hi_mom0.beam.jtok(hi_freq) / u.Jy).to(u.K * u.km / u.s) * np.cos(gal.inclination)
hi_mom0_5beam_Kkms = (hi_mom0_5beam * hi_mom0.beam.jtok(hi_freq) / u.Jy).to(u.K * u.km / u.s) * np.cos(gal.inclination)

hi_surfdens_2beam = hi_mass_conversion * hi_mom0_2beam_Kkms
hi_surfdens_5beam = hi_mass_conversion * hi_mom0_5beam_Kkms

# Compare distributions of the HI at different resolutions:

twocolumn_twopanel_figure()

plt.subplot(121)

orig = astro_hist(hi_surfdens.value[hi_mask], bins='knuth', alpha=0.5,
                  label='80 pc')

twobeam = astro_hist(hi_surfdens_2beam.value[hi_mask & np.isfinite(hi_surfdens_2beam)],
                     bins=orig[1], alpha=0.5,
                     label='160 pc')

fivebeam = astro_hist(hi_surfdens_5beam.value[hi_mask & np.isfinite(hi_surfdens_5beam)],
                      bins=orig[1], alpha=0.5,
                      label='380 pc')

plt.legend(loc='best')
plt.xlabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ / pc$^2$)")

plt.grid()

# Make CDFs
plt.subplot(122)

plt.plot((orig[1][1:] + orig[1][:-1]) * 0.5,
         np.cumsum(orig[0]) / np.sum(orig[0]),
         drawstyle='steps-mid')
plt.plot((orig[1][1:] + orig[1][:-1]) * 0.5,
         np.cumsum(twobeam[0]) / np.sum(twobeam[0]),
         drawstyle='steps-mid')
plt.plot((orig[1][1:] + orig[1][:-1]) * 0.5,
         np.cumsum(fivebeam[0]) / np.sum(fivebeam[0]),
         drawstyle='steps-mid')

plt.xlabel(r"$\Sigma_{\rm HI}$ (M$_{\odot}$ / pc$^2$)")

plt.grid()

plt.tight_layout()

plt.savefig(osjoin(fig_path, "hi_surfdens_resolution_change_hist.pdf"))
plt.savefig(osjoin(fig_path, "hi_surfdens_resolution_change_hist.png"))
plt.close()

# Now recompute the pressures with the lower res HI

pext_on_gmc_2beam = \
    Pext_sum(co_surfdens,
             hi_surfdens_2beam,
             stellar_surfdens_rep,
             H_mol=100 * u.pc,
             H_star=300 * u.pc,
             sigma_HI=7 * u.km / u.s)

pext_on_gmc_5beam = \
    Pext_sum(co_surfdens,
             hi_surfdens_5beam,
             stellar_surfdens_rep,
             H_mol=100 * u.pc,
             H_star=300 * u.pc,
             sigma_HI=7 * u.km / u.s)

# How does the HI pressure change?

twocolumn_twopanel_figure()

ax1 = plt.subplot(121)

_ = astro_hist(np.log10(pext_on_gmc[3][co_mask].value),
               bins='knuth', alpha=0.5,
               label='80 pc')
_ = astro_hist(np.log10(pext_on_gmc_2beam[3][co_mask].value),
               bins='knuth', alpha=0.5,
               label='160 pc')
_ = astro_hist(np.log10(pext_on_gmc_5beam[3][co_mask].value),
               bins='knuth', alpha=0.5,
               label='380 pc')

plt.xlabel(r"log Pressure (K cm$^{-3}$)")
plt.grid()
plt.text(0.05, 0.85, "HI on GMC", transform=ax1.transAxes,
         bbox={"boxstyle": "square", "facecolor": "w",
               "edgecolor": "gray"})
plt.legend(loc='lower left')

ax2 = plt.subplot(122)

_ = astro_hist(np.log10(pext_on_gmc[4][co_mask].value),
               bins='knuth', alpha=0.5,
               label='80 pc')
_ = astro_hist(np.log10(pext_on_gmc_2beam[4][co_mask].value),
               bins='knuth', alpha=0.5,
               label='160 pc')
_ = astro_hist(np.log10(pext_on_gmc_5beam[4][co_mask].value),
               bins='knuth', alpha=0.5,
               label='380 pc')

plt.grid()
plt.text(0.05, 0.85, "Stellar on HI", transform=ax2.transAxes,
         bbox={"boxstyle": "square", "facecolor": "w",
               "edgecolor": "gray"})

plt.xlabel(r"log Pressure (K cm$^{-3}$)")

plt.tight_layout()

plt.savefig(osjoin(fig_path, "pext_hiterms_resolution_comparison.pdf"))
plt.savefig(osjoin(fig_path, "pext_hiterms_resolution_comparison.png"))
plt.close()

# Load in tables from M33-ish mass galaxies in Jiayi's sample

phangs_data_dir = os.path.expanduser("~/bigdata/ekoch/PHANGS/")

tab_ngc4540 = Table.read(osjoin(phangs_data_dir, "NGC4540_phys_1000pc.fits"))
tab_ngc5042 = Table.read(osjoin(phangs_data_dir, "NGC5042_phys_1000pc.fits"))
tab_ngc5068 = Table.read(osjoin(phangs_data_dir, "NGC5068_phys_1000pc.fits"))

# Compare total pressure to M33. The other galaxies are avg'd over
# kpc scales (weighted by Gauss kernel; Leroy+16)
# Compare the CDFs

return_finite = lambda tab, y: tab[y][np.isfinite(tab[y])]

twocolumn_figure()

fig, ax = plt.subplots(2, 2, sharey=True, sharex=True)

# Total pressure

vals_m33, bins_m33 = \
    histogram(np.log10(pext_on_gmc_5beam[0][co_mask].value),
              bins='knuth')

ax[0, 0].plot((bins_m33[1:] + bins_m33[:-1]) * 0.5,
              np.cumsum(vals_m33) / vals_m33.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='M33 (380 pc HI)')

p_4540 = return_finite(tab_ngc4540, "<P_DE_pix_90pc>")
vals_4540, bins_4540 = \
    histogram(np.log10(p_4540),
              bins='knuth')
ax[0, 0].plot((bins_4540[1:] + bins_4540[:-1]) * 0.5,
              np.cumsum(vals_4540) / vals_4540.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 4540')

p_5042 = return_finite(tab_ngc5042, "<P_DE_pix_90pc>")
vals_5042, bins_5042 = \
    histogram(np.log10(p_5042),
              bins='knuth')
ax[0, 0].plot((bins_5042[1:] + bins_5042[:-1]) * 0.5,
              np.cumsum(vals_5042) / vals_5042.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5042')

p_5068 = return_finite(tab_ngc5068, "<P_DE_pix_90pc>")
vals_5068, bins_5068 = \
    histogram(np.log10(p_5068),
              bins='knuth')
ax[0, 0].plot((bins_5068[1:] + bins_5068[:-1]) * 0.5,
              np.cumsum(vals_5068) / vals_5068.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5068')

ax[0, 0].text(0.05, 0.85, "Total",
              transform=ax[0, 0].transAxes,
              bbox={"boxstyle": "square", "facecolor": "w",
                    "edgecolor": "gray"})
ax[0, 0].grid()
ax[0, 0].legend(loc='lower left')

# GMC pressure

vals_m33, bins_m33 = \
    histogram(np.log10(pext_on_gmc_5beam[1][co_mask].value),
              bins='knuth')

ax[1, 0].plot((bins_m33[1:] + bins_m33[:-1]) * 0.5,
              np.cumsum(vals_m33) / vals_m33.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='M33 (380 pc HI)')

p_4540 = return_finite(tab_ngc4540, "<P_DE^molg_pix_90pc>")
vals_4540, bins_4540 = \
    histogram(np.log10(p_4540),
              bins='knuth')
ax[1, 0].plot((bins_4540[1:] + bins_4540[:-1]) * 0.5,
              np.cumsum(vals_4540) / vals_4540.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 4540')

p_5042 = return_finite(tab_ngc5042, "<P_DE^molg_pix_90pc>")
vals_5042, bins_5042 = \
    histogram(np.log10(p_5042),
              bins='knuth')
ax[1, 0].plot((bins_5042[1:] + bins_5042[:-1]) * 0.5,
              np.cumsum(vals_5042) / vals_5042.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5042')

p_5068 = return_finite(tab_ngc5068, "<P_DE^molg_pix_90pc>")
vals_5068, bins_5068 = \
    histogram(np.log10(p_5068),
              bins='knuth')
ax[1, 0].plot((bins_5068[1:] + bins_5068[:-1]) * 0.5,
              np.cumsum(vals_5068) / vals_5068.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5068')

ax[1, 0].text(0.05, 0.85, "GMC Press.",
              transform=ax[1, 0].transAxes,
              bbox={"boxstyle": "square", "facecolor": "w",
                    "edgecolor": "gray"})
ax[1, 0].grid()

# Stellar pressure on GMC

vals_m33, bins_m33 = \
    histogram(np.log10(pext_on_gmc_5beam[2][co_mask].value),
              bins='knuth')

ax[0, 1].plot((bins_m33[1:] + bins_m33[:-1]) * 0.5,
              np.cumsum(vals_m33) / vals_m33.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='M33 (380 pc HI)')

p_4540 = return_finite(tab_ngc4540, "<P_DE^starg_pix_90pc>")
vals_4540, bins_4540 = \
    histogram(np.log10(p_4540),
              bins='knuth')
ax[0, 1].plot((bins_4540[1:] + bins_4540[:-1]) * 0.5,
              np.cumsum(vals_4540) / vals_4540.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 4540')

p_5042 = return_finite(tab_ngc5042, "<P_DE^starg_pix_90pc>")
vals_5042, bins_5042 = \
    histogram(np.log10(p_5042),
              bins='knuth')
ax[0, 1].plot((bins_5042[1:] + bins_5042[:-1]) * 0.5,
              np.cumsum(vals_5042) / vals_5042.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5042')

p_5068 = return_finite(tab_ngc5068, "<P_DE^starg_pix_90pc>")
vals_5068, bins_5068 = \
    histogram(np.log10(p_5068),
              bins='knuth')
ax[0, 1].plot((bins_5068[1:] + bins_5068[:-1]) * 0.5,
              np.cumsum(vals_5068) / vals_5068.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5068')

ax[0, 1].text(0.05, 0.85, "Stellar on GMC",
              transform=ax[0, 1].transAxes,
              bbox={"boxstyle": "square", "facecolor": "w",
                    "edgecolor": "gray"})
ax[0, 1].grid()

# Stellar pressure on HI + HI self-grav

vals_m33, bins_m33 = \
    histogram(np.log10(pext_on_gmc_5beam[3][co_mask].value +
                       pext_on_gmc_5beam[4][co_mask].value),
              bins='knuth')

ax[1, 1].plot((bins_m33[1:] + bins_m33[:-1]) * 0.5,
              np.cumsum(vals_m33) / vals_m33.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='M33 (380 pc HI)')

p_4540 = return_finite(tab_ngc4540, "P_DE^atom")
vals_4540, bins_4540 = \
    histogram(np.log10(p_4540),
              bins='knuth')
ax[1, 1].plot((bins_4540[1:] + bins_4540[:-1]) * 0.5,
              np.cumsum(vals_4540) / vals_4540.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 4540')

p_5042 = return_finite(tab_ngc5042, "P_DE^atom")
vals_5042, bins_5042 = \
    histogram(np.log10(p_5042),
              bins='knuth')
ax[1, 1].plot((bins_5042[1:] + bins_5042[:-1]) * 0.5,
              np.cumsum(vals_5042) / vals_5042.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5042')

p_5068 = return_finite(tab_ngc5068, "P_DE^atom")
vals_5068, bins_5068 = \
    histogram(np.log10(p_5068),
              bins='knuth')
ax[1, 1].plot((bins_5068[1:] + bins_5068[:-1]) * 0.5,
              np.cumsum(vals_5068) / vals_5068.sum(),
              drawstyle='steps-mid',
              alpha=0.5,
              label='NGC 5068')

ax[1, 1].text(0.05, 0.85, "HI Press.",
              transform=ax[1, 1].transAxes,
              bbox={"boxstyle": "square", "facecolor": "w",
                    "edgecolor": "gray"})
ax[1, 1].grid()


fig.text(0.5, 0.04, r"log Pressure (K cm$^{-3}$)",
         ha='center', va='center')
fig.text(0.06, 0.5, r'CDF',
         ha='center', va='center', rotation='vertical')

plt.savefig(osjoin(fig_path, "pext_m33_phangs_cdf_comparison.pdf"))
plt.savefig(osjoin(fig_path, "pext_m33_phangs_cdf_comparison.png"))
plt.close()

