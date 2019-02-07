from galaxies import Galaxy
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
from astropy.table import Table
from scipy import interpolate
import astropy.units as u
from astropy.io import fits

g = Galaxy('M33')
t = Table.read('rad.out.csv', format='ascii.csv')
R = t['r'].data * 0.0122 * 1e3
vrot = t['Vt'].data
R = np.linspace(0, 1.2e4, 1e5)
r_pc = R
H = 100 * u.pc
rn = R / (1.2e4)
n = 0.56
vrot = 110.0 * rn / (0.33333 + 0.666666 * rn**n)**(3 / (2 * n))
t, c, k = interpolate.splrep(r_pc, vrot, s=0, k=3)
knots = 8
# Cubic interpolation of vrot(R).
vrot_sm = interpolate.BSpline(t, c, k, extrapolate=True)
Nsteps = 10000
R = np.linspace(R.min(), R.max(), Nsteps)

# plt.plot(R, vrot_sm(R))
# plt.plot(r_pc, vrot)

dVdR = interpolate.splev(R, (t, c, k), der=1)
dOmdlnR = (dVdR - vrot_sm(R) / R) * u.km / u.s / u.pc
Bdens = 0.6 * (8.0e-6)**2 / (8 * np.pi) * u.erg / u.cm**3
edot_mri = (-Bdens * dOmdlnR).to(u.erg / u.cm**3 / u.s)

sfr = fits.getdata('m33_sfr.fuv24.fits') * np.cos(g.inclination)
sfr_hdr = fits.getheader('m33_sfr.fuv24.fits')

rgal = g.radius(header=sfr_hdr).to(u.pc)
stat, edges, count = ss.binned_statistic(rgal.value.ravel(),
                                         sfr.ravel(),
                                         statistic=np.nanmedian,
                                         bins=np.linspace(0, 12000, 50))
centers = 0.5 * (edges[0:-1] + edges[1:])
snr = stat / u.kpc**2 / u.yr * 1e-2  # SNR per area
eps_SNR = 0.1  # Efficience of supernova
edot_snr = (snr * 1e51 * u.erg * eps_SNR / H).to(u.erg / u.cm**3 / u.s)
edot_obs_min = 1.5 * (8 * u.M_sun / u.pc**2
                      * (6 * u.km / u.s)**3
                      / H**2).to(u.erg / u.cm**3 / u.s)
edot_obs_max = 1.5 * (8 * u.M_sun / u.pc**2
                      * (12 * u.km / u.s)**3
                      / H**2).to(u.erg / u.cm**3 / u.s)

fig, ax = plt.subplots(1, 1)
ax.semilogy(centers / 1e3, edot_snr, label='SN', color='red')
ax.semilogy(R / 1e3, edot_mri, label='MRI', color='blue', linestyle='--')
ones = np.ones_like(R)
ax.fill_between(R / 1e3, ones * edot_obs_min.value,
                ones * edot_obs_max.value,
                label='Obs.', color='gray', alpha=0.35)
ax.set_xlabel(r'$R_{\mathrm{gal}}\ (\mathrm{kpc})$')
ax.set_ylabel(r'$\dot{u}\ (\mathrm{erg\ s^{-1}\ cm^{-3}})$')
ax.set_ylim(1e-28, 1e-24)
ax.set_xlim(0, 8)
ax.legend(loc='upper right')
fig.set_size_inches(4, 3)
plt.tight_layout()
fig.savefig('m33_energy.pdf')
