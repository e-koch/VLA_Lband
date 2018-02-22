
'''
Examine the line shapes that result from a mixture of thermal line width
Gaussian components.

I'm getting the moments of the normal mixture model from:
http://onlyvix.blogspot.ca/2011/05/mixture-of-normal-formulas-for-skew-and.html
which uses the derivation from:
GENERATING DAILY CHANGES IN MARKET VARIABLES
USING A MULTIVARIATE MIXTURE OF NORMAL DISTRIBUTIONS
Jin Wang 2001.

'''

import numpy as np
import astropy.units as u
import astropy.constants as co
import matplotlib.pyplot as plt
from scipy.special import erf

from cube_analysis.spectral_stacking_models import fit_hwhm


def mix_mu(mu_samps, f=None):

    if f is not None:
        assert mu_samps.size == f.size
        assert f.sum() == 1.
    else:
        f = np.ones_like(mu_samps) / float(mu_samps.size)

    return np.sum(f * mu_samps)


def mix_variance(mu_samps, var_samps, f=None, mu=None):

    if mu is None:
        mu = mix_mu(mu_samps, f)

    if f is not None:
        assert mu_samps.size == f.size
        assert f.sum() == 1.
    else:
        f = np.ones_like(mu_samps) / float(mu_samps.size)

    term1 = np.sum(f * (var_samps + mu_samps**2)) - mu**2

    return term1 - mu**2


def mix_skew(mu_samps, var_samps, f=None, mu=None, var=None):

    if mu is None:
        mu = mix_mu(mu_samps, f)

    if var is None:
        var = mix_variance(mu_samps, var_samps, f, mu=mu)

    if f is not None:
        assert mu_samps.size == f.size
        assert f.sum() == 1.
    else:
        f = np.ones_like(mu_samps) / float(mu_samps.size)

    term1 = np.sum(f * (mu_samps - mu) * (3 * var + (mu_samps - mu)**2))

    return term1 / var**1.5


def mix_kurt(mu_samps, var_samps, f=None, mu=None, var=None):

    if mu is None:
        mu = mix_mu(mu_samps, f)

    if var is None:
        var = mix_variance(mu_samps, var_samps, f, mu=mu)

    terma = 3 * var_samps**2
    termb = 6 * (mu_samps - mu)**2 * var_samps
    termc = (mu_samps - mu)**4

    if f is not None:
        assert mu_samps.size == f.size
        assert f.sum() == 1.
    else:
        f = np.ones_like(mu_samps) / float(mu_samps.size)

    term1 = np.sum(f * (terma + termb + termc))

    return term1 / var**2


def lwidth(T):
    return np.sqrt((co.k_B * T) / (1.4 * co.m_p)).to(u.km / u.s)


def gauss(mu, std):
    return lambda x: np.exp(-(x - mu)**2 / (2 * std**2))


def pearson7_profile(mu, std, kurt):
    '''
    A Pearson VII distribution with kurtosis as a free parameter.
    https://en.wikipedia.org/wiki/Pearson_distribution
    '''

    excess_kurt = kurt - 3.

    m = 2.5 + 3 / excess_kurt

    alpha = std * np.sqrt(2 * m - 3)

    return lambda x: (1 + ((x - mu) / alpha)**2) ** (-m)


def create_profile(vels, mu_samps, std_samps):
    '''
    Create the equivalent stacked profile shape.
    '''

    yvals = np.zeros_like(vels)

    for mu, std in zip(mu_samps, std_samps):
        yvals += gauss(mu, std)(vels)

    return yvals / yvals.max()


# Draw a few thousand samples according to different distributions for the mean
# and sigma.
nsamp = 10000

vels = np.linspace(-20, 20, 40 / 0.2)
# Uniform temperature from 100 to 1e4 K

peak_centers = np.random.normal(0, 0.1, nsamp)
rotation_centers = np.random.normal(0, 2, nsamp)

uniform_temps = np.random.uniform(100., 10000., nsamp) * u.K
# uniform_temps = np.random.uniform(1000., 10000., nsamp) * u.K

uniform_lwidths = lwidth(uniform_temps).value

# Estimate the stacked properties.
stack_lwidth = np.sqrt(mix_variance(peak_centers, uniform_lwidths**2))
stack_lwidth_rot = np.sqrt(mix_variance(rotation_centers, uniform_lwidths**2))

uniform_profile = create_profile(vels, peak_centers, uniform_lwidths)
uniform_profile_rot = create_profile(vels, rotation_centers, uniform_lwidths)

# What is the equivalent width we would measure with the HWHM technique?
params, _, labels, hwhm_model = fit_hwhm(vels, uniform_profile)
params_rot, _, _, hwhm_model_rot = fit_hwhm(vels, uniform_profile_rot)

# Calculate skewness and kurtosis
uniform_skew = mix_skew(peak_centers, uniform_lwidths**2)
uniform_skew_rot = mix_skew(rotation_centers, uniform_lwidths**2)

uniform_kurt = mix_kurt(peak_centers, uniform_lwidths**2)
uniform_kurt_rot = mix_kurt(rotation_centers, uniform_lwidths**2)

plt.subplot(211)
plt.plot(vels, uniform_profile, label='Mixture')
plt.plot(vels, gauss(0, stack_lwidth)(vels), label='Reference')
plt.plot(vels, hwhm_model(vels), label='HWHM Model')
plt.legend(frameon=True)
plt.grid()

plt.subplot(212)
plt.plot(vels, uniform_profile_rot, label='Mixture')
plt.plot(vels, gauss(0, stack_lwidth_rot)(vels), label='Reference')
plt.plot(vels, hwhm_model_rot(vels), label='HWHM Model')
plt.grid()

print("Peak Velocity Params")
print("sigma_mix = {0}".format(stack_lwidth))
print(params)
print("Rotation Velocity Params")
print("sigma_mix = {0}".format(stack_lwidth_rot))
print(params_rot)

print(argh)

# Now try injecting a more realistic temperature distribution.
# Using f_CNM ~ 0.2, f_unstab ~ 0.2 and f_WNM ~0.6 from Murray et al. 2015

f_CNM = 0.2
f_unstab = 0.2
f_WNM = 0.6

T_CNM = np.random.uniform(40, 500, int(f_CNM * nsamp))
T_unstab = np.random.uniform(500, 5000, int(f_unstab * nsamp))
T_WNM = np.random.uniform(5000, 10000, int(f_WNM * nsamp))

phase_temps = np.hstack([T_CNM, T_unstab, T_WNM])

phase_lwidths = lwidth(phase_temps * u.K).value

# Estimate the stacked properties.
stack_lwidth_phase = np.sqrt(mix_variance(peak_centers, phase_lwidths**2))
stack_lwidth_phase_rot = np.sqrt(mix_variance(rotation_centers, phase_lwidths**2))

phase_profile = create_profile(vels, peak_centers, phase_lwidths)
phase_profile_rot = create_profile(vels, rotation_centers, phase_lwidths)

# What is the equivalent width we would measure with the HWHM technique?
params_phase = fit_hwhm(vels, phase_profile)[0]
params_phase_rot = fit_hwhm(vels, phase_profile_rot)[0]

# Calculate skewness and kurtosis
phase_skew = mix_skew(peak_centers, phase_lwidths**2)
phase_skew_rot = mix_skew(rotation_centers, phase_lwidths**2)

phase_kurt = mix_kurt(peak_centers, phase_lwidths**2)
phase_kurt_rot = mix_kurt(rotation_centers, phase_lwidths**2)

plt.subplot(211)
plt.plot(vels, phase_profile, label='Mixture')
plt.plot(vels, gauss(0, stack_lwidth_phase)(vels), label='Reference')
plt.plot(vels, gauss(0, sigma_hwhm_phase)(vels), label='HWHM Model')
plt.legend(frameon=True)
plt.grid()

plt.subplot(212)
plt.plot(vels, phase_profile_rot, label='Mixture')
plt.plot(vels, gauss(0, stack_lwidth_phase_rot)(vels), label='Reference')
plt.plot(vels, gauss(0, sigma_hwhm_phase_rot)(vels), label='HWHM Model')
plt.grid()

print("Peak Velocity Params")
print("sigma_mix = {0}; sigma_hwhm = {1}; f_wings = {2}; mix_skew = {3}"
      "; mix_kurt = {4}".format(stack_lwidth_phase, sigma_hwhm_phase,
                                f_wings_phase, phase_skew, phase_kurt))
print("Rotation Velocity Params")
print("sigma_mix = {0}; sigma_hwhm = {1}; f_wings = {2}; mix_skew = {3}"
      "; mix_kurt = {4}".format(stack_lwidth_phase_rot, sigma_hwhm_phase_rot,
                                f_wings_rot,
                                phase_skew_rot, phase_kurt_rot))

# What about the canonical temps?

# T_CNM_canon = np.abs(np.random.normal(200, 50, int(0.2 * nsamp)))
T_CNM_canon = np.random.uniform(40, 200, int(0.2 * nsamp))
# T_WNM_canon = np.abs(np.random.normal(8000, 50, int(0.8 * nsamp)))
T_WNM_canon = np.random.uniform(5000, 1000, int(0.8 * nsamp))

canon_temps = np.hstack([T_CNM_canon, T_WNM_canon])

canon_lwidths = lwidth(canon_temps * u.K).value

# Estimate the stacked properties.
stack_lwidth_canon = np.sqrt(mix_variance(peak_centers, canon_lwidths**2))
stack_lwidth_canon_rot = np.sqrt(mix_variance(rotation_centers, canon_lwidths**2))

canon_profile = create_profile(vels, peak_centers, canon_lwidths)
canon_profile_rot = create_profile(vels, rotation_centers, canon_lwidths)

# What is the equivalent width we would measure with the HWHM technique?
sigma_hwhm_canon, f_wings_canon = fit_hwhm(vels, canon_profile)[0]
sigma_hwhm_canon_rot, f_wings_canon_rot = fit_hwhm(vels, canon_profile_rot)[0]

# Calculate skewness and kurtosis
canon_skew = mix_skew(peak_centers, canon_lwidths**2)
canon_skew_rot = mix_skew(rotation_centers, canon_lwidths**2)

canon_kurt = mix_kurt(peak_centers, canon_lwidths**2)
canon_kurt_rot = mix_kurt(rotation_centers, canon_lwidths**2)

plt.subplot(211)
plt.plot(vels, canon_profile, label='Mixture')
plt.plot(vels, gauss(0, stack_lwidth_canon)(vels), label='Reference')
plt.plot(vels, gauss(0, sigma_hwhm_canon)(vels), label='HWHM Model')
plt.legend(frameon=True)
plt.grid()

plt.subplot(212)
plt.plot(vels, canon_profile_rot, label='Mixture')
plt.plot(vels, gauss(0, stack_lwidth_canon_rot)(vels), label='Reference')
plt.plot(vels, gauss(0, sigma_hwhm_canon_rot)(vels), label='HWHM Model')
plt.grid()

print("Peak Velocity Params")
print("sigma_mix = {0}; sigma_hwhm = {1}; f_wings = {2}; mix_skew = {3}"
      "; mix_kurt = {4}".format(stack_lwidth_canon, sigma_hwhm_canon,
                                f_wings_canon, canon_skew, canon_kurt))
print("Rotation Velocity Params")
print("sigma_mix = {0}; sigma_hwhm = {1}; f_wings = {2}; mix_skew = {3}"
      "; mix_kurt = {4}".format(stack_lwidth_canon_rot, sigma_hwhm_canon_rot,
                                f_wings_rot,
                                canon_skew_rot, canon_kurt_rot))


# Compare the uniform and the canonical temperature distributions.

plt.plot(vels, uniform_profile, label='Uniform T')
plt.plot(vels, phase_profile, label='Phase T')
plt.plot(vels, canon_profile, label='Canon T')
plt.legend(frameon=True)
plt.grid()


# Finally, what is the effect of multi-component spectra?
# Let's assume a small fraction (5%) have multi-components
# Of these, they will have an average offset that scatters
# with a std of the typical line width

f_multi = 0.5

peak_cent_single = np.random.normal(0, 0.1, int((1 - f_multi) * nsamp))

# To sample the multi-components, let's just assume they will increase the
# dispersion of the peaks, and so the dominant component adds preferentially
# to the line wings
# peak_cent_multi = np.random.normal(0, 6, int(f_multi * nsamp))
peak_cent_multi = np.random.uniform(-6, 6, int(f_multi * nsamp))

peak_centers_wmulti = np.hstack([peak_cent_single, peak_cent_multi])

# Estimate the stacked properties.
stack_lwidth_multi = np.sqrt(mix_variance(peak_centers_wmulti,
                                          uniform_lwidths**2))

uniform_profile_multi = create_profile(vels, peak_centers_wmulti,
                                       uniform_lwidths)

# What is the equivalent width we would measure with the HWHM technique?
sigma_hwhm_multi, f_wings_multi = fit_hwhm(vels, uniform_profile_multi)[0]

# Calculate skewness and kurtosis
uniform_skew_multi = mix_skew(peak_centers_wmulti, uniform_lwidths**2)

uniform_kurt_multi = mix_kurt(peak_centers_wmulti, uniform_lwidths**2)

plt.plot(vels, uniform_profile, label='Uniform T')
plt.plot(vels, uniform_profile_multi, label='Uniform T w/ multi')
plt.legend(frameon=True)
plt.grid()
