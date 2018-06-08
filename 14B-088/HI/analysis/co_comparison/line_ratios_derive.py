
'''
Compare HI/CO line width ratio expected from a turbulent medium.
'''

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import astropy.units as u
import astropy.constants as c

from paths import allfigs_path

from plotting_styles import (default_figure, twocolumn_twopanel_figure,
                             onecolumn_figure)


default_figure()


def lwidth(T, mu):
    return np.sqrt(c.k_B * T / (mu * c.m_p)).to(u.km / u.s)


def size_lwidth_turb(l, l_s, T, mu, index=0.5):
    '''
    Assume shock-dominated, so index is 1/2.
    '''

    c_s = lwidth(T, mu)

    return ((l / l_s)**index * c_s).to(u.km / u.s)


def sonic_length(sigma, l, T, mu, index=0.5):
    '''
    Sonic length assuming sigma is the total line width (turb + therm).
    Turbulent s-l index still 1/2.
    '''

    c_s = lwidth(T, mu)

    return l / ((sigma / c_s)**2 - 1)**(2 * index)


def size_lwidth(l, l_s, T, mu, index=0.5):
    '''
    Full size-line width relation.
    '''

    c_s = lwidth(T, mu)

    turb_sigma = size_lwidth_turb(l, l_s, T, mu, index=index)

    return np.sqrt(turb_sigma**2 + c_s**2)


# From stacking
sigma_HI = 6.6 * u.km / u.s
sigma_CO = 4.5 * u.km / u.s

# ASSUME that the turbulent field is the same through the cloud and the
# envelope.
# Set the size-line width relation based on sigma_CO since it's dominated by
# turbulence

T_CO = 10 * u.K
sigma_turb = np.sqrt(sigma_CO**2 - lwidth(T_CO, 2.3)**2)

# Set T_HI assuming this holds for both
c_HI = np.sqrt(sigma_HI**2 - sigma_turb**2)
T_HI = (c_HI**2 * 1.4 * c.m_p / c.k_B).to(u.K)
# T_HI = 5000 * u.K
# About 4000 K. Reasonable.

scales = np.arange(0.05, 1500, 0.1) * u.pc

beam_size = 80 * u.pc

l_HI = sonic_length(sigma_HI, beam_size, T_HI, 1.4)
l_CO = sonic_length(sigma_CO, beam_size, T_CO, 2.3)

sigma_COs = size_lwidth(scales, l_CO, T_CO, 2.3)
sigma_HIs = size_lwidth(scales, l_HI, T_HI, 1.4)

# What if the HI index is closer to 1/3?
l_HI_one3 = sonic_length(sigma_HI, beam_size, T_HI, 1.4, index=1 / 3.)
sigma_HIs_one3 = size_lwidth(scales, l_HI_one3, T_HI, 1.4, index=1 / 3.)

twocolumn_twopanel_figure()
plt.subplot(121)
plt.loglog(scales.value, sigma_COs.value, label='CO')
plt.loglog(scales.value, sigma_HIs.value, label='HI')
plt.loglog(scales.value, sigma_HIs_one3.value, label='HI - 1/3 index')
plt.axhline(0.19 * (100 / l_CO.value)**0.5)
plt.errorbar([80, 160, 380], [4.5, 6.0, 7.25], yerr=[1.3] * 3,
             fmt='D', label='CO Width')
plt.errorbar([80, 160, 380], [6.6, 8.0, 8.9], yerr=[0.1] * 3,
             fmt='o', label='HI Width')
plt.grid()
plt.legend(frameon=True)
plt.xlabel("Scale (pc)")
plt.ylabel(r"$\sigma$ (km/s)")
plt.subplot(122)
plt.loglog(scales.value, sigma_COs.value / sigma_HIs.value,
           label='HI Index = 1/2')
plt.loglog(scales.value, sigma_COs.value / sigma_HIs_one3.value,
           label='HI Index = 1/3')
plt.errorbar([80, 160, 380], [4.5 / 6.6, 6.0 / 8.0, 7.25 / 8.9],
             yerr=[0.2, 0.16, 0.15],
             fmt='D')
plt.plot([40], [0.23], 'D')
plt.axhline(1., color='k', linestyle='--')
plt.grid()
plt.legend(frameon=True)
plt.xlabel("Scale (pc)")
plt.ylabel("Ratio")

plt.tight_layout()

# plt.savefig(allfigs_path("co_vs_hi/size_linewidth_peaksub_ratio.png"))
# plt.savefig(allfigs_path("co_vs_hi/size_linewidth_peaksub_ratio.pdf"))
plt.close()

onecolumn_figure()

# Long legend requires a slightly larger figure
plt.figure(figsize=(4.5, 3.24))

# Need some extra colors
sb.set_palette('colorblind', n_colors=8)

cpal = sb.color_palette()

plt.loglog(scales.value, sigma_COs.value, color=cpal[0],
           label='CO - 1/2 index')
plt.loglog(scales.value, sigma_HIs.value, "--", color=cpal[1],
           label='HI - 1/2 index')
plt.loglog(scales.value, sigma_HIs_one3.value, "-.", color=cpal[1],
           label='HI - 1/3 index')
plt.errorbar([80, 160, 380], [4.5, 6.0, 7.25], yerr=[1.3] * 3, color=cpal[2],
             fmt='D', label='M33 CO(2-1) Width')
plt.errorbar([80, 160, 380], [6.6, 8.0, 8.9], yerr=[0.1] * 3, color=cpal[2],
             fmt='o', label='M33 HI Width')
# Show the mean value from Fukui+09
plt.errorbar([40], [14.1 / 2.35], yerr=[3.3 / 2.35], fmt='s',
             label='LMC HI Avg.',
             color=cpal[3])
plt.errorbar([40], [4.6 / 2.35], yerr=[1.6 / 2.35], fmt='X',
             label='LMC CO(1-0) Avg.',
             color=cpal[3])

# Narrow line width in M31 from Caldu-Primo + Schruba 2016
plt.errorbar([85], [7.5 / 2.35], yerr=[0.4], fmt='X', label='M31 CO(1-0) Stack',
             color=cpal[5])

# NGC 2403 CO and HI from Caldu-Primo+2013
ngc2403_widths = np.array([[12.3, 15.3, 0.8],
                           [12.0, 16.4, 0.7],
                           [12.4, 14.2, 0.9],
                           [12.4, 11.1, 1.1],
                           [12.2, 9.9, 1.2],
                           [12.2, 10.3, 1.2],
                           [12.1, 9.1, 1.3],
                           [11.1, 8.8, 1.3],
                           [10.7, 9.1, 1.2],
                           [11.2, 9.4, 1.2],
                           [10.6, 8.2, 1.3],
                           [10.7, 7.6, 1.4],
                           [10.5, 8.4, 1.2],
                           [10.8, 8.7, 1.2],
                           [10.6, 8.3, 1.3],
                           [10.3, 6.7, 1.5],
                           [10.8, 9.6, 1.1],
                           [10.2, 9.8, 1.0],
                           [10.8, 7.6, 1.4],
                           [10.8, 7.2, 1.5],
                           [9.7, 12.2, 0.8]])
hi_err = 1.0
co_err = 3.0

# Exclude inner points that are likely affected by beam smearing
plt.errorbar([200], ngc2403_widths[3:, 0].mean(),
             yerr=np.sqrt(hi_err**2 + ngc2403_widths[3:, 0].std()**2),
             fmt='s', label='NGC 2403 HI Avg.',
             color=cpal[7])
plt.errorbar([200], ngc2403_widths[3:, 1].mean(),
             yerr=np.sqrt(co_err**2 + ngc2403_widths[3:, 1].std()**2),
             fmt='X', label='NGC 2403 CO(2-1) Avg.',
             color=cpal[7])

print("NGC 2403 HI: {0}+/-{1}".format(ngc2403_widths[3:, 0].mean(),
                                      np.sqrt(hi_err**2 +
                                              ngc2403_widths[3:, 0].std()**2)))
print("NGC 2403 CO: {0}+/-{1}".format(ngc2403_widths[3:, 1].mean(),
                                      np.sqrt(co_err**2 +
                                              ngc2403_widths[3:, 1].std()**2)))

plt.grid()
plt.legend(frameon=True)
plt.xlabel("Scale (pc)")
plt.ylabel(r"$\sigma$ (km/s)")
plt.tight_layout()

plt.savefig(allfigs_path("co_vs_hi/size_linewidth_peaksub_comparison.png"))
plt.savefig(allfigs_path("co_vs_hi/size_linewidth_peaksub_comparison.pdf"))
plt.close()

# sigma_CO_cutoff = np.zeros_like(sigma_COs)
# sigma_CO_cutoff[scales.value < 100] = size_lwidth(scales, l_CO, T_CO, 2.3)[scales.value < 100]
# sigma_CO_cutoff[scales.value >= 100] = lwidth(10 * u.K, 2.3) * (100 / l_CO.value)**(0.5)

# # p.plot(scales.value, np.sqrt(sigma_COs.value**2 + (0.02 * scales.value)**2))
# plt.loglog(scales.value, sigma_COs)
# plt.loglog(scales.value, sigma_CO_cutoff)

# plt.loglog(scales, sigma_CO_cutoff / sigma_HIs)
