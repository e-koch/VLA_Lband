
'''
Compare HI/CO line width ratio expected from a turbulent medium.
'''

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import astropy.units as u
import astropy.constants as c


from plotting_styles import default_figure


default_figure()


def lwidth(T, mu):
    return np.sqrt(c.k_B * T / (mu * c.m_p)).to(u.km / u.s)


def size_lwidth_turb(l, l_s, T, mu):
    '''
    Assume shock-dominated, so index is 1/2.
    '''

    c_s = lwidth(T, mu)

    return (np.sqrt(l / l_s) * c_s).to(u.km / u.s)


def sonic_length(sigma, l, T, mu):
    '''
    Sonic length assuming sigma is the total line width (turb + therm).
    Turbulent s-l index still 1/2.
    '''

    c_s = lwidth(T, mu)

    return l / ((sigma / c_s)**2 - 1)


def size_lwidth(l, l_s, T, mu):
    '''
    Full size-line width relation.
    '''

    c_s = lwidth(T, mu)

    turb_sigma = size_lwidth_turb(l, l_s, T, mu)

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
# About 4000 K. Reasonable.

scales = np.arange(0.05, 1500, 0.1) * u.pc

beam_size = 80 * u.pc

l_HI = sonic_length(sigma_HI, beam_size, T_HI, 1.4)
l_CO = sonic_length(sigma_CO, beam_size, T_CO, 2.3)

sigma_COs = size_lwidth(scales, l_CO, T_CO, 2.3)
sigma_HIs = size_lwidth(scales, l_HI, T_HI, 1.4)

plt.loglog(scales.value, sigma_COs.value, label='CO')
plt.loglog(scales.value, sigma_HIs.value, label='HI')
