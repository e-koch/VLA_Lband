
'''
Compare different estimates of the line width and how finite sampling
affects the line shape.
'''

import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from astropy.modeling import models, fitting
from astropy.convolution import convolve_fft
from scipy.special import erf
from scipy.interpolate import InterpolatedUnivariateSpline
from os.path import join as osjoin
from scipy.optimize import curve_fit

from plotting_styles import onecolumn_figure, twocolumn_figure
from paths import allfigs_path

fig_folder = allfigs_path("finite_sampled_fitting")


def gauss_model(vels, amp=1., stddev=10., mean=0.):
    return models.Gaussian1D(amplitude=amp, stddev=stddev, mean=mean)(vels)


def gauss_weighted_avg(a, b, amp=1., stddev=10., mean=0.):

    SQRT2 = np.sqrt(2)

    erf_term = erf((mean - a) / (SQRT2 * stddev)) - \
        erf((mean - b) / (SQRT2 * stddev))

    return (amp * stddev / (b - a)) * erf_term * np.sqrt(np.pi / 2.)


def gauss_model_chanweight(vels, amp=1., stddev=10., mean=0.):

    # Assume the channels are equally spaced.
    chan_width = np.abs(np.diff(vels[:2])[0])

    low_vels = vels - chan_width / 2.
    high_vels = vels + chan_width / 2.

    vals = np.zeros_like(vels)

    for i, (low, high) in enumerate(zip(low_vels, high_vels)):
        vals[i] = gauss_weighted_avg(low, high, amp, stddev, mean)

    return vals


# Empirical measures of width


def equiv_gauss_width(integrated, peak):
    return integrated / (np.sqrt(2 * np.pi) * peak)


def moment2_width(vels, spec):

    cent = np.sum(spec * vels) / np.sum(spec)

    return np.sqrt(np.sum(spec * (vels - cent)**2) / np.sum(spec))


def find_hwhm(vels, spectrum, interp_factor=10):
    '''
    Return the equivalent Gaussian sigma based on the HWHM positions.
    '''

    halfmax = spectrum.max() * 0.5

    # Model the spectrum with a spline
    # x values must be increasing for the spline, so flip if needed.
    spec_for_interp, vels_for_interp = reorder_spectra(vels, spectrum)

    interp1 = InterpolatedUnivariateSpline(vels_for_interp,
                                           spec_for_interp - halfmax, k=3)

    fwhm_points = interp1.roots()
    if len(fwhm_points) < 2:
        raise ValueError("Found less than 2 roots!")
    # Only keep the min/max if there are multiple
    fwhm_points = (min(fwhm_points), max(fwhm_points))

    fwhm = fwhm_points[1] - fwhm_points[0]

    # Convert to equivalent Gaussian sigma
    sigma = fwhm / np.sqrt(8 * np.log(2))

    # Upsample in velocity to estimate the peak position
    interp_factor = float(interp_factor)
    chan_size = np.diff(vels_for_interp[:2])[0] / interp_factor
    upsamp_vels = np.linspace(vels_for_interp.min(),
                              vels_for_interp.max() + 0.9 * chan_size,
                              vels_for_interp.size * interp_factor)
    upsamp_spec = interp1(upsamp_vels)
    peak_velocity = upsamp_vels[np.argmax(upsamp_spec)]

    return sigma, fwhm_points, peak_velocity


def reorder_spectra(vels, spectrum):
    spec_for_interp = spectrum if np.diff(vels[:2])[0] > 0 else spectrum[::-1]
    vels_for_interp = vels if np.diff(vels[:2])[0] > 0 else vels[::-1]

    return spec_for_interp, vels_for_interp


def leroy16_deconv(chan_width, r=0.26):
    k = 0.0 + 0.47 * r - 0.23 * r**2 - 0.16 * r**3 + 0.43 * r**4
    return chan_width * (1.0 + 1.18 * k + 10.4 * k**2) / np.sqrt(2 * np.pi)


def hanning_kernel(r=0.26):
    '''
    Correlation for M33 from Sun+18
    '''
    k = 0.0 + 0.47 * r - 0.23 * r**2 - 0.16 * r**3 + 0.43 * r**4
    return np.array([k, 1 - 2 * k, k])


def convolve_gauss_model_chanweight(vels, kernel, amp=1., stddev=10., mean=0.):
    return convolve_fft(gauss_model_chanweight(vels, amp=amp, stddev=stddev,
                                               mean=mean), kernel)


def fit_discrete_model(vels, spec, p0=(1., 10., 0.)):
    '''
    Use curve_fit to fit the broadened model.
    '''
    return curve_fit(lambda vels, amp, stddev, mean:
                     gauss_model_chanweight(vels, amp=amp,
                                            stddev=stddev,
                                            mean=mean),
                     vels, spec, p0=p0)


def fit_convolve_model(vels, spec, kernel, p0=(1., 10., 0.)):
    '''
    Use curve_fit to fit the broadened model.
    '''
    return curve_fit(lambda vels, amp, stddev, mean:
                     convolve_gauss_model_chanweight(vels, kernel, amp=amp,
                                                     stddev=stddev,
                                                     mean=mean),
                     vels, spec, p0=p0)


# Plot an example Gaussian, samples at the bin centre, finite samples in
# 1 sigma channels, and the effect of broadening.

sigma = 10.
r = 0.26

m33_kernel = hanning_kernel(r=r)

vels_5 = np.append(- np.arange(0., 61., 5.)[::-1][:-1], np.arange(0., 61., 5.))
vels_10 = np.append(- np.arange(0., 61., 10.)[::-1][:-1],
                    np.arange(0., 61., 10.))
vels_01 = np.append(- np.arange(0., 60., 0.1)[::-1][:-1],
                    np.arange(0., 60., 0.1))

plt.plot(vels_01 / sigma, gauss_model(vels_01), label='Model',
         drawstyle='steps-mid', linewidth=4, alpha=0.4)
plt.plot(vels_10 / sigma, gauss_model(vels_10), label='Model at \nbin centre',
         drawstyle='steps-mid')
plt.plot(vels_10 / sigma, gauss_model_chanweight(vels_10), '--',
         label='Finite-sampling', drawstyle='steps-mid')
plt.plot(vels_10 / sigma, convolve_gauss_model_chanweight(vels_10, m33_kernel),
         '-.', label='With response \nfunction', drawstyle='steps-mid')
plt.legend(frameon=True, loc='upper right')
plt.grid()
plt.xlim([-7.3, 7.3])
plt.ylabel("Amplitude")
plt.xlabel(r"$N \, \sigma$")
plt.tight_layout()

plt.savefig(osjoin(fig_folder, "example_model_chan_1sigma.png"))
plt.savefig(osjoin(fig_folder, "example_model_chan_1sigma.pdf"))
plt.close()

plt.plot((vels_01 - 5.) / sigma, gauss_model(vels_01 - 5.), label='Model',
         drawstyle='steps-mid', linewidth=4, alpha=0.4)
plt.plot((vels_10 - 5.) / sigma, gauss_model(vels_10 - 5.),
         label='Model at \nbin centre',
         drawstyle='steps-mid')
plt.plot((vels_10 - 5.) / sigma, gauss_model_chanweight(vels_10 - 5.), '--',
         label='Finite-sampling', drawstyle='steps-mid')
plt.plot((vels_10 - 5.) / sigma,
         convolve_gauss_model_chanweight(vels_10 - 5., m33_kernel),
         '-.', label='With response \nfunction', drawstyle='steps-mid')
plt.legend(frameon=True, loc='upper right')
plt.grid()
plt.xlim([-7.3, 7.3])
plt.ylabel("Amplitude")
plt.xlabel(r"$N \, \sigma$")
plt.tight_layout()

plt.savefig(osjoin(fig_folder, "example_model_chan_1sigma_chanoffset.png"))
plt.savefig(osjoin(fig_folder, "example_model_chan_1sigma_chanoffset.pdf"))
plt.close()

# Now compare the different width estimates as a function of channel width

chan_widths = np.array([0.1, 0.2, 0.5, 1., 2.5, 5., 10., 15., 20., 30.])

amp = 1.
stddev = 10.
mean = 0.

init_mod = models.Gaussian1D(amplitude=amp, stddev=stddev, mean=mean)

stddev_actual = []
stddev_weight = []
stddev_weight_deconv = []
stddev_weight_formod = []
stddev_weight_equiv = []
stddev_weight_equiv_deconv = []
stddev_weight_mom2 = []
stddev_weight_mom2_deconv = []
stddev_weight_hwhm = []
stddev_weight_hwhm_deconv = []

stddev_weight_resp = []
stddev_weight_resp_formod = []
stddev_weight_resp_deconv = []
stddev_weight_resp_equiv = []
stddev_weight_resp_equiv_deconv = []
stddev_weight_resp_mom2 = []
stddev_weight_resp_mom2_deconv = []
stddev_weight_resp_hwhm = []
stddev_weight_resp_hwhm_deconv = []

for width in chan_widths:

    # Best case: The centre is falls at the middle of a channel
    vels = np.append(- np.arange(0., 60. + width / 2., width)[::-1][:-1],
                     np.arange(0., 60. + width / 2., width))

    gauss_spec = gauss_model(vels)
    gauss_weight_spec = gauss_model_chanweight(vels)

    fitter = fitting.LevMarLSQFitter()

    gauss_fit = fitter(init_mod, vels, gauss_spec)
    stddev_actual.append(gauss_fit.stddev)

    gauss_weight_fit = fitter(init_mod, vels, gauss_weight_spec)
    stddev_weight.append(gauss_weight_fit.stddev)

    stddev_weight_deconv.append(np.sqrt(gauss_weight_fit.stddev**2 -
                                        leroy16_deconv(width, r=0.)**2))

    # Correct fit
    gauss_weight_disc_fit = fit_discrete_model(vels, gauss_weight_spec)
    stddev_weight_formod.append(gauss_weight_disc_fit[0][1])

    # Now the empirical width measures
    equiv_width = equiv_gauss_width(gauss_weight_spec.sum() * width, gauss_weight_spec.max())

    stddev_weight_equiv.append(equiv_width)
    stddev_weight_equiv_deconv.append(np.sqrt(equiv_width**2 - leroy16_deconv(width, r=0.)**2))

    mom2_width = moment2_width(vels, gauss_weight_spec)

    stddev_weight_mom2.append(mom2_width)
    stddev_weight_mom2_deconv.append(np.sqrt(mom2_width**2 - leroy16_deconv(width, r=0.)**2))

    # This fails if the channel width is beyond the FWHM. Expected, but the function does not
    # gracefully fail
    if width < stddev * 2.35:
        hwhm_width = find_hwhm(vels, gauss_weight_spec)[0]
    else:
        hwhm_width = np.NaN

    stddev_weight_hwhm.append(hwhm_width)
    stddev_weight_hwhm_deconv.append(np.sqrt(hwhm_width**2 - leroy16_deconv(width, r=0.)**2))

    # And now for the response-convolved spectrum
    gauss_weight_spec_conv = convolve_gauss_model_chanweight(vels, m33_kernel)

    # Correct fit
    gauss_weight_resp_fit = fit_convolve_model(vels, gauss_weight_spec_conv,
                                               m33_kernel)
    stddev_weight_resp_formod.append(gauss_weight_resp_fit[0][1])

    # Normal fit
    gauss_weight_fit = fitter(init_mod, vels, gauss_weight_spec_conv)
    stddev_weight_resp.append(gauss_weight_fit.stddev)

    stddev_weight_resp_deconv.append(np.sqrt(gauss_weight_fit.stddev**2 -
                                             leroy16_deconv(width, r=r)**2))

    # Now the empirical width measures
    equiv_width = equiv_gauss_width(gauss_weight_spec_conv.sum() * width, gauss_weight_spec_conv.max())

    stddev_weight_resp_equiv.append(equiv_width)
    stddev_weight_resp_equiv_deconv.append(np.sqrt(equiv_width**2 -
                                                   leroy16_deconv(width, r=0.26)**2))

    mom2_width = moment2_width(vels, gauss_weight_spec_conv)

    stddev_weight_resp_mom2.append(mom2_width)
    stddev_weight_resp_mom2_deconv.append(np.sqrt(mom2_width**2 -
                                                  leroy16_deconv(width, r=0.26)**2))

    # This fails if the channel width is beyond the FWHM. Expected, but the function does not
    # gracefully fail
    if width < stddev * 2.35:
        hwhm_width = find_hwhm(vels, gauss_weight_spec_conv)[0]
    else:
        hwhm_width = np.NaN

    stddev_weight_resp_hwhm.append(hwhm_width)
    stddev_weight_resp_hwhm_deconv.append(np.sqrt(hwhm_width**2 -
                                                  leroy16_deconv(width, r=0.26)**2))


# Now run everything with the line centre set at 1/2 channel

init_mod = models.Gaussian1D(amplitude=amp, stddev=stddev, mean=mean)

stddev_actual_offset = []
stddev_weight_offset = []
stddev_weight_deconv_offset = []
stddev_weight_formod_offset = []
stddev_weight_equiv_offset = []
stddev_weight_equiv_deconv_offset = []
stddev_weight_mom2_offset = []
stddev_weight_mom2_deconv_offset = []
stddev_weight_hwhm_offset = []
stddev_weight_hwhm_deconv_offset = []

stddev_weight_resp_offset = []
stddev_weight_resp_formod_offset = []
stddev_weight_resp_deconv_offset = []
stddev_weight_resp_equiv_offset = []
stddev_weight_resp_equiv_deconv_offset = []
stddev_weight_resp_mom2_offset = []
stddev_weight_resp_mom2_deconv_offset = []
stddev_weight_resp_hwhm_offset = []
stddev_weight_resp_hwhm_deconv_offset = []

for width in chan_widths:

    # Best case: The centre is falls at the middle of a channel
    vels = np.append(- np.arange(0., 60. + width / 2., width)[::-1][:-1],
                     np.arange(0., 60. + width / 2., width))

    gauss_spec = gauss_model(vels, mean=width / 2.)
    gauss_weight_spec = gauss_model_chanweight(vels, mean=width / 2.)

    fitter = fitting.LevMarLSQFitter()

    gauss_fit = fitter(init_mod, vels, gauss_spec)
    stddev_actual_offset.append(gauss_fit.stddev)

    gauss_weight_fit = fitter(init_mod, vels, gauss_weight_spec)
    stddev_weight_offset.append(gauss_weight_fit.stddev)

    stddev_weight_deconv_offset.append(np.sqrt(gauss_weight_fit.stddev**2 -
                                               leroy16_deconv(width, r=0.)**2))

    # Correct fit
    gauss_weight_disc_fit = fit_discrete_model(vels, gauss_weight_spec)
    stddev_weight_formod_offset.append(gauss_weight_disc_fit[0][1])

    # Now the empirical width measures
    equiv_width = equiv_gauss_width(gauss_weight_spec.sum() * width, gauss_weight_spec.max())

    stddev_weight_equiv_offset.append(equiv_width)
    stddev_weight_equiv_deconv_offset.append(np.sqrt(equiv_width**2 - leroy16_deconv(width, r=0.)**2))

    mom2_width = moment2_width(vels, gauss_weight_spec)

    stddev_weight_mom2_offset.append(mom2_width)
    stddev_weight_mom2_deconv_offset.append(np.sqrt(mom2_width**2 - leroy16_deconv(width, r=0.)**2))

    # This fails if the channel width is beyond the FWHM. Expected, but the function does not
    # gracefully fail
    if width < stddev * 2.35:
        hwhm_width = find_hwhm(vels, gauss_weight_spec)[0]
    else:
        hwhm_width = np.NaN

    stddev_weight_hwhm_offset.append(hwhm_width)
    stddev_weight_hwhm_deconv_offset.append(np.sqrt(hwhm_width**2 -
                                             leroy16_deconv(width, r=0.)**2))

    # And now for the response-convolved spectrum
    gauss_weight_spec_conv = \
        convolve_gauss_model_chanweight(vels, m33_kernel, mean=width / 2.)

    # Correct fit
    gauss_weight_resp_fit = fit_convolve_model(vels, gauss_weight_spec_conv,
                                               m33_kernel)
    stddev_weight_resp_formod_offset.append(gauss_weight_resp_fit[0][1])

    # Normal fit
    gauss_weight_fit = fitter(init_mod, vels, gauss_weight_spec_conv)
    stddev_weight_resp_offset.append(gauss_weight_fit.stddev)

    stddev_weight_resp_deconv_offset.append(np.sqrt(gauss_weight_fit.stddev**2 -
                                             leroy16_deconv(width, r=r)**2))

    # Now the empirical width measures
    equiv_width = equiv_gauss_width(gauss_weight_spec_conv.sum() * width, gauss_weight_spec_conv.max())

    stddev_weight_resp_equiv_offset.append(equiv_width)
    stddev_weight_resp_equiv_deconv_offset.append(np.sqrt(equiv_width**2 -
                                                   leroy16_deconv(width, r=0.26)**2))

    mom2_width = moment2_width(vels, gauss_weight_spec_conv)

    stddev_weight_resp_mom2_offset.append(mom2_width)
    stddev_weight_resp_mom2_deconv_offset.append(np.sqrt(mom2_width**2 -
                                                  leroy16_deconv(width, r=0.26)**2))

    # This fails if the channel width is beyond the FWHM. Expected, but the function does not
    # gracefully fail
    if width < stddev * 2.35:
        hwhm_width = find_hwhm(vels, gauss_weight_spec_conv)[0]
    else:
        hwhm_width = np.NaN

    stddev_weight_resp_hwhm_offset.append(hwhm_width)
    stddev_weight_resp_hwhm_deconv_offset.append(np.sqrt(hwhm_width**2 -
                                                  leroy16_deconv(width, r=0.26)**2))

# Arrays are nice
stddev_actual = np.array(stddev_actual)
stddev_weight = np.array(stddev_weight)
stddev_weight_deconv = np.array(stddev_weight_deconv)
stddev_weight_formod = np.array(stddev_weight_formod)
stddev_weight_equiv = np.array(stddev_weight_equiv)
stddev_weight_equiv_deconv = np.array(stddev_weight_equiv_deconv)
stddev_weight_mom2 = np.array(stddev_weight_mom2)
stddev_weight_mom2_deconv = np.array(stddev_weight_mom2_deconv)
stddev_weight_hwhm = np.array(stddev_weight_hwhm)
stddev_weight_hwhm_deconv = np.array(stddev_weight_hwhm_deconv)

stddev_weight_resp = np.array(stddev_weight_resp)
stddev_weight_resp_formod = np.array(stddev_weight_resp_formod)
stddev_weight_resp_deconv = np.array(stddev_weight_resp_deconv)
stddev_weight_resp_equiv = np.array(stddev_weight_resp_equiv)
stddev_weight_resp_equiv_deconv = np.array(stddev_weight_resp_equiv_deconv)
stddev_weight_resp_mom2 = np.array(stddev_weight_resp_mom2)
stddev_weight_resp_mom2_deconv = np.array(stddev_weight_resp_mom2_deconv)
stddev_weight_resp_hwhm = np.array(stddev_weight_resp_hwhm)
stddev_weight_resp_hwhm_deconv = np.array(stddev_weight_resp_hwhm_deconv)

stddev_actual_offset = np.array(stddev_actual_offset)
stddev_weight_offset = np.array(stddev_weight_offset)
stddev_weight_deconv_offset = np.array(stddev_weight_deconv_offset)
stddev_weight_formod_offset = np.array(stddev_weight_formod_offset)
stddev_weight_equiv_offset = np.array(stddev_weight_equiv_offset)
stddev_weight_equiv_deconv_offset = np.array(stddev_weight_equiv_deconv_offset)
stddev_weight_mom2_offset = np.array(stddev_weight_mom2_offset)
stddev_weight_mom2_deconv_offset = np.array(stddev_weight_mom2_deconv_offset)
stddev_weight_hwhm_offset = np.array(stddev_weight_hwhm_offset)
stddev_weight_hwhm_deconv_offset = np.array(stddev_weight_hwhm_deconv_offset)

stddev_weight_resp_offset = np.array(stddev_weight_resp_offset)
stddev_weight_resp_formod_offset = np.array(stddev_weight_resp_formod_offset)
stddev_weight_resp_deconv_offset = np.array(stddev_weight_resp_deconv_offset)
stddev_weight_resp_equiv_offset = np.array(stddev_weight_resp_equiv_offset)
stddev_weight_resp_equiv_deconv_offset = np.array(stddev_weight_resp_equiv_deconv_offset)
stddev_weight_resp_mom2_offset = np.array(stddev_weight_resp_mom2_offset)
stddev_weight_resp_mom2_deconv_offset = np.array(stddev_weight_resp_mom2_deconv_offset)
stddev_weight_resp_hwhm_offset = np.array(stddev_weight_resp_hwhm_offset)
stddev_weight_resp_hwhm_deconv_offset = np.array(stddev_weight_resp_hwhm_deconv_offset)

# Summarize the various tests into one 4-panel figure
twocolumn_figure()

cpal = sb.color_palette()

fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)

# First panel is channel centered w/ finite sampling

# axs[0, 0].axvline(1., color=cpal[0])
axs[0, 0].plot(sigma / chan_widths, stddev_weight / sigma, color=cpal[0],
               label='Fit to Gaussian')
axs[0, 0].plot(sigma / chan_widths, stddev_weight_deconv / sigma, '--', color=cpal[0])
axs[0, 0].plot(sigma / chan_widths, stddev_weight_equiv / sigma, color=cpal[1],
               label='Equivalent Width')
axs[0, 0].plot(sigma / chan_widths, stddev_weight_equiv_deconv / sigma, '--',
               color=cpal[1])
axs[0, 0].plot(sigma / chan_widths, stddev_weight_mom2 / sigma, color=cpal[2],
               label='Moment 2')
axs[0, 0].plot(sigma / chan_widths, stddev_weight_mom2_deconv / sigma, '--',
               color=cpal[2])
axs[0, 0].plot(sigma / chan_widths, stddev_weight_hwhm / sigma, color=cpal[3],
               label='HWHM')
axs[0, 0].plot(sigma / chan_widths, stddev_weight_hwhm_deconv / sigma, '--',
               color=cpal[3])
axs[0, 0].plot(sigma / chan_widths, stddev_weight_formod / sigma, color=cpal[5],
               label='Fit to Discrete Model')
axs[0, 0].legend(frameon=True, fontsize=9)
axs[0, 0].grid()
axs[0, 0].set_title("Peak at channel centre")
# axs[0, 0].text(0.5, 2.5, "Peak at channel centre",
#               bbox={"boxstyle": "square", "facecolor": "w"})

# First column, second row is the convolved model
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp / sigma, color=cpal[0],
               label='Fit to Gaussian')
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_deconv / sigma, '--', color=cpal[0])
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_equiv / sigma, color=cpal[1],
               label='Equivalent Width')
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_equiv_deconv / sigma, '--',
               color=cpal[1])
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_mom2 / sigma, color=cpal[2],
               label='Moment 2')
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_mom2_deconv / sigma, '--',
               color=cpal[2])
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_hwhm / sigma, color=cpal[3],
               label='HWHM')
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_hwhm_deconv / sigma, '--',
               color=cpal[3])
axs[1, 0].plot(sigma / chan_widths, stddev_weight_resp_formod / sigma, color=cpal[5],
               label='Fit to Discrete Model')
axs[1, 0].grid()

# Second row are with half a channel offset
axs[0, 1].plot(sigma / chan_widths, stddev_weight_offset / sigma, color=cpal[0],
               label='Fit to Gaussian')
axs[0, 1].plot(sigma / chan_widths, stddev_weight_deconv_offset / sigma, '--', color=cpal[0])
axs[0, 1].plot(sigma / chan_widths, stddev_weight_equiv_offset / sigma, color=cpal[1],
               label='Equivalent Width')
axs[0, 1].plot(sigma / chan_widths, stddev_weight_equiv_deconv_offset / sigma, '--',
               color=cpal[1])
axs[0, 1].plot(sigma / chan_widths, stddev_weight_mom2_offset / sigma, color=cpal[2],
               label='Moment 2')
axs[0, 1].plot(sigma / chan_widths, stddev_weight_mom2_deconv_offset / sigma, '--',
               color=cpal[2])
axs[0, 1].plot(sigma / chan_widths, stddev_weight_hwhm_offset / sigma, color=cpal[3],
               label='HWHM')
axs[0, 1].plot(sigma / chan_widths, stddev_weight_hwhm_deconv_offset / sigma, '--',
               color=cpal[3])
axs[0, 1].plot(sigma / chan_widths, stddev_weight_formod_offset / sigma, color=cpal[5],
               label='Fit to Discrete Model')
axs[0, 1].grid()
axs[0, 1].set_title("Peak at channel edge")

# First column, second row is the convolved model
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_offset / sigma, color=cpal[0],
               label='Fit to Gaussian')
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_deconv_offset / sigma, '--', color=cpal[0])
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_equiv_offset / sigma, color=cpal[1],
               label='Equivalent Width')
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_equiv_deconv_offset / sigma, '--',
               color=cpal[1])
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_mom2_offset / sigma, color=cpal[2],
               label='Moment 2')
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_mom2_deconv_offset / sigma, '--',
               color=cpal[2])
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_hwhm_offset / sigma, color=cpal[3],
               label='HWHM')
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_hwhm_deconv_offset / sigma, '--',
               color=cpal[3])
axs[1, 1].plot(sigma / chan_widths, stddev_weight_resp_formod_offset / sigma, color=cpal[5],
               label='Fit to Discrete Model')
axs[1, 1].grid()

axs[1, 1].set_xlim([0.2, 3])
# axs[1, 1].set_ylim([8, 15])
# axs[1].set_xlabel(r"Channels per $\sigma$")

fig.text(0.5, 0.02, r'Channels per $\sigma$', ha='center', va='center')
fig.text(0.03, 0.5, r'Measured Width / $\sigma$',
         ha='center', va='center', rotation='vertical')

fig.text(0.95, 0.7, 'Ideal\nResponse Function', ha='center',
         va='center', rotation='vertical')
fig.text(0.95, 0.3, 'Hanning\nResponse Function', ha='center',
         va='center', rotation='vertical')

fig.subplots_adjust(hspace=0.04, wspace=0.03)

fig.savefig(osjoin(fig_folder, "width_recovery_comparison.png"))
fig.savefig(osjoin(fig_folder, "width_recovery_comparison.pdf"))

plt.close()
