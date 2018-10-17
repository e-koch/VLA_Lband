
'''
Take the cloud positions from the Corbelli+20176 catalogue and calculate
the two-point correlation function. This demonstrates that there is excess
correlation up to scales of ~2 kpc (the scale of the CO emission), which
results in a non-zero slope on large-scale, as found in Combes+2012.
'''

from astropy.table import Table
import astropy.units as u
from astroML.correlation import bootstrap_two_point
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from os.path import join as osjoin

from paths import data_path, allfigs_path
from constants import ang_to_phys
from plotting_styles import onecolumn_figure, onecolumn_twopanel_figure

tab = Table.read(osjoin(data_path, "Corbelli_17_catalogues/J_A+A_601_A146_table5.dat.fits"))

# These are in deg to match the catalogue
# Choose the smallest bin to be larger than most cloud radii (~120 pc).
# Scales below this have large negative correlations.
bins = np.linspace(30. / 3600., 700. / 3600., 20)

data = np.vstack([tab['RAdeg'], tab['DEdeg']])

phys_bins = ang_to_phys(0.5 * (bins[1:] + bins[:-1]) * u.deg).to(u.kpc)

test, test_err = bootstrap_two_point(data.T, bins, Nbootstrap=2000)

onecolumn_figure()

plt.errorbar(phys_bins.value, test, yerr=test_err,
             drawstyle='steps-mid')
plt.axhline(0., linestyle='--', color='k')
plt.grid()
plt.xlabel("Separation (kpc)")
plt.ylabel("Correlation")

plt.tight_layout()

plt.savefig(allfigs_path("CO21_properties/corbelli_cloudloc_2ptcorr.pdf"))
plt.savefig(allfigs_path("CO21_properties/corbelli_cloudloc_2ptcorr.png"))
plt.close()

# Run one more with the Landy-Szalay estimator against a random sample

# radeg_R = np.random.uniform(low=tab['RAdeg'].min(), high=tab['RAdeg'].max())
# decdeg_R = np.random.uniform(low=tab['DEdeg'].min(), high=tab['DEdeg'].max())

# data_R = np.vstack([radeg_R, decdeg_R])

test_ls, test_err_ls = bootstrap_two_point(data.T, bins, Nbootstrap=2000,
                                           method='landy-szalay')


onecolumn_twopanel_figure()

fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)

ax[0].errorbar(phys_bins.value, test, yerr=test_err, label='Standard',
               drawstyle='steps-mid')
ax[0].axhline(0., linestyle='--', color='k')
ax[0].legend(frameon=True)
ax[0].grid()

ax[1].errorbar(phys_bins.value, test_ls, yerr=test_err_ls,
               label='Landy-Szalay',
               drawstyle='steps-mid')
ax[1].axhline(0., linestyle='--', color='k')
ax[1].legend(frameon=True)
ax[1].grid()
ax[1].set_xlabel("Separation (kpc)")

fig.savefig(allfigs_path("CO21_properties/corbelli_cloudloc_2ptcorr_wlandyszalay.pdf"))
fig.savefig(allfigs_path("CO21_properties/corbelli_cloudloc_2ptcorr_wlandyszalay.png"))
plt.close()
