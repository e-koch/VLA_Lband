
'''
How well do the GBT and VLA data match in the overlap regions?
'''

from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import os
from radio_beam import Beam
import astropy.units as u
import numpy as np
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.base.model import GenericLikelihoodModel
from astropy.visualization import hist

from cube_analysis.feather_cubes import feather_compare_cube

from paths import fourteenB_HI_data_path, data_path, allfigs_path
from constants import hi_freq
from plotting_styles import onecolumn_figure, default_figure, twocolumn_figure


vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))

gbt_path = os.path.join(data_path, "GBT")
gbt_registered_cube = SpectralCube.read(os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid_registered.fits"))

beam_fwhm = lambda diam: ((1.18 * hi_freq.to(u.cm, u.spectral())) / diam.to(u.cm)) * u.rad

# The shortest baseline in the 14B-088 data is ~44 m.
las = (hi_freq.to(u.cm, u.spectral()) / (44 * u.m)).to(u.arcsec, u.dimensionless_angles())

# Return the distributions of ratios for all channels.
radii, ratios, high_pts, low_pts = \
    feather_compare_cube(vla_cube, gbt_registered_cube, las,
                         num_cores=6,
                         lowresfwhm=beam_fwhm(90 * u.m))

radii_85, ratios_85, high_pts_85, low_pts_85 = \
    feather_compare_cube(vla_cube, gbt_registered_cube, las,
                         num_cores=6,
                         lowresfwhm=beam_fwhm(85 * u.m))

radii_875, ratios_875, high_pts_875, low_pts_875 = \
    feather_compare_cube(vla_cube, gbt_registered_cube, las,
                         num_cores=6,
                         lowresfwhm=beam_fwhm(87.5 * u.m))

radii_95, ratios_95, high_pts_95, low_pts_95 = \
    feather_compare_cube(vla_cube, gbt_registered_cube, las,
                         num_cores=6,
                         lowresfwhm=beam_fwhm(95 * u.m))

# Test #1 -- what are the slopes between k and ratio per channel?


def sentheil_perchan(xvals, yvals, alpha=0.85):

    slope = np.empty((len(xvals)))
    upper_uncert = np.empty((len(xvals)))
    lower_uncert = np.empty((len(xvals)))

    for i, (xval, yval) in enumerate(zip(xvals, yvals)):

        out = stats.theilslopes(yval, x=xval, alpha=alpha)

        slope[i] = out[0]
        upper_uncert[i] = out[3] - out[0]
        lower_uncert[i] = out[0] - out[2]

    return slope, lower_uncert, upper_uncert


slope, lower_uncert, upper_uncert = \
    sentheil_perchan([rad.to(u.arcmin)**2 for rad in radii], ratios)
slope_85, lower_uncert_85, upper_uncert_85 = \
    sentheil_perchan([rad.to(u.arcmin)**2 for rad in radii_85], ratios_85)
slope_875, lower_uncert_875, upper_uncert_875 = \
    sentheil_perchan([rad.to(u.arcmin)**2 for rad in radii_875], ratios_875)
slope_95, lower_uncert_95, upper_uncert_95 = \
    sentheil_perchan([rad.to(u.arcmin)**2 for rad in radii_95], ratios_95)

chans = np.arange(len(radii))

# For visualizing, do a LOWESS smoothing
slope_lowess = lowess(slope, chans, frac=0.1, is_sorted=True,
                      return_sorted=False)
slope_lowess_85 = lowess(slope_85, chans, frac=0.1, is_sorted=True,
                         return_sorted=False)
slope_lowess_875 = lowess(slope_875, chans, frac=0.1, is_sorted=True,
                          return_sorted=False)
slope_lowess_95 = lowess(slope_95, chans, frac=0.1, is_sorted=True,
                         return_sorted=False)

twocolumn_figure()

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
ax[0, 0].errorbar(chans, slope_85,
                  yerr=[lower_uncert_85, upper_uncert_85], label="85 m",
                  alpha=0.5)
ax[0, 0].plot(chans, slope_lowess_85)
ax[0, 0].axhline(0, linestyle='--')
ax[0, 0].text(0, 0.015, "{}".format(np.round(beam_fwhm(85 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 0].set_ylabel(r"Slope (arcmin$^{-2}$)")
ax[0, 0].grid(True)

ax[0, 1].errorbar(chans, slope_875,
                  yerr=[lower_uncert_875, upper_uncert_875], label="87.5 m",
                  alpha=0.5)
ax[0, 1].plot(chans, slope_lowess_875)
ax[0, 1].text(0, 0.015, "{}".format(np.round(beam_fwhm(87.5 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 1].axhline(0, linestyle='--')
ax[0, 1].grid(True)

ax[1, 0].errorbar(chans, slope, yerr=[lower_uncert, upper_uncert],
                  label="90 m", alpha=0.5)
ax[1, 0].plot(chans, slope_lowess)
ax[1, 0].axhline(0, linestyle='--')
ax[1, 0].text(0, 0.015, "{}".format(np.round(beam_fwhm(90 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 0].set_ylabel(r"Slope (arcmin$^{-2}$)")
ax[1, 0].set_xlabel("Channel")
ax[1, 0].grid(True)

ax[1, 1].errorbar(chans, slope_95,
                  yerr=[lower_uncert_95, upper_uncert_95], label="95 m",
                  alpha=0.5)
ax[1, 1].plot(chans, slope_lowess_95)
ax[1, 1].axhline(0, linestyle='--')
ax[1, 1].text(0, 0.015, "{}".format(np.round(beam_fwhm(95 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 1].set_xlabel("Channel")
ax[1, 1].grid(True)

fig.tight_layout()

fig.savefig(allfigs_path("angres_vs_ratios_vla_gbt.png"))
fig.savefig(allfigs_path("angres_vs_ratios_vla_gbt.pdf"))

plt.close()

# Make a separate plots for the 87.5 m (9.8') beam

onecolumn_figure()
plt.errorbar(chans, slope_875,
             yerr=[lower_uncert_875, upper_uncert_875], label="87.5 m",
             alpha=0.5)
plt.plot(chans, slope_lowess_875)
plt.axhline(0, linestyle='--')
plt.grid(True)
plt.ylabel(r"Slope (arcmin$^{-2}$)")
plt.xlabel("Channel")
plt.tight_layout()
plt.savefig(allfigs_path("angres_vs_ratios_vla_gbt_9.8arcmin.png"))
plt.savefig(allfigs_path("angres_vs_ratios_vla_gbt_9.8arcmin.pdf"))
plt.close()

# Test #2 -- Sen-Theil fit to the low and high res points to get scale factor

# Do per channel
sf_slope, sf_slope_llim, sf_slope_hlim = sentheil_perchan(low_pts, high_pts)
sf_slope_85, sf_slope_llim_85, sf_slope_hlim_85 = \
    sentheil_perchan(low_pts_85, high_pts_85)
sf_slope_875, sf_slope_llim_875, sf_slope_hlim_875 = \
    sentheil_perchan(low_pts_875, high_pts_875)
sf_slope_95, sf_slope_llim_95, sf_slope_hlim_95 = \
    sentheil_perchan(low_pts_95, high_pts_95)

twocolumn_figure()

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
ax[0, 0].errorbar(chans, sf_slope_85,
                  yerr=[sf_slope_llim_85, sf_slope_hlim_85], label="85 m",
                  alpha=0.5)
# ax[0, 0].plot(chans, slope_lowess_85)
ax[0, 0].axhline(1, linestyle='--')
ax[0, 0].text(0, 1.2, "{}".format(np.round(beam_fwhm(85 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 0].set_ylabel(r"Scale Factor")
ax[0, 0].grid(True)

ax[0, 1].errorbar(chans, sf_slope_875,
                  yerr=[sf_slope_llim_875, sf_slope_hlim_875], label="87.5 m",
                  alpha=0.5)
# ax[0, 1].plot(chans, slope_lowess_875)
ax[0, 1].text(0, 1.2, "{}".format(np.round(beam_fwhm(87.5 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 1].axhline(1, linestyle='--')
ax[0, 1].grid(True)

ax[1, 0].errorbar(chans, sf_slope, yerr=[sf_slope_llim, sf_slope_hlim],
                  label="90 m", alpha=0.5)
# ax[1, 0].plot(chans, slope_lowess)
ax[1, 0].axhline(1, linestyle='--')
ax[1, 0].text(0, 1.2, "{}".format(np.round(beam_fwhm(90 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 0].set_ylabel(r"Scale Factor")
ax[1, 0].set_xlabel("Channel")
ax[1, 0].grid(True)

ax[1, 1].errorbar(chans, sf_slope_95,
                  yerr=[sf_slope_llim_95, sf_slope_hlim_95], label="95 m",
                  alpha=0.5)
# ax[1, 1].plot(chans, slope_lowess_95)
ax[1, 1].axhline(1, linestyle='--')
ax[1, 1].text(0, 1.2, "{}".format(np.round(beam_fwhm(95 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 1].set_xlabel("Channel")
ax[1, 1].grid(True)

fig.tight_layout()

fig.savefig(allfigs_path("scalefactors_fitted_vla_gbt.png"))
fig.savefig(allfigs_path("scalefactors_fitted_vla_gbt.pdf"))

plt.close()

# Test #3 -- Cauchy fit to the ratios to get the scale factor


class Likelihood(GenericLikelihoodModel):

    # Get the number of parameters from shapes.
    # Add one for scales, since we're assuming loc is frozen.
    # Keeping loc=0 is appropriate for log-normal models.
    nparams = 1 if stats.cauchy.shapes is None else \
        len(stats.cauchy.shapes.split(",")) + 1

    def loglike(self, params):
        if np.isnan(params).any():
            return - np.inf

        loglikes = \
            stats.cauchy.logpdf(self.endog, *params[:-2],
                                scale=params[-1],
                                loc=params[-2])
        if not np.isfinite(loglikes).all():
            return - np.inf
        else:
            return loglikes.sum()


def cauchy_fitter(ratios, verbose=False):

    init_fit = stats.cauchy.fit(np.log(ratios))

    mle_model = Likelihood(np.log(ratios))
    fitted_model = mle_model.fit(init_fit, method='nm', disp=0)
    fitted_model.df_model = len(ratios)
    fitted_model.df_resid = len(ratios) - 2

    if verbose:
        _ = hist(np.log(ratios), bins='scott', normed=True, alpha=0.5)
        xvals = np.arange(-3, 3, 0.01)
        plt.plot(xvals, stats.cauchy.pdf(xvals, *fitted_model.params))

    return fitted_model.params[0], fitted_model.bse[0], fitted_model


# Fit a Cauchy distribution to the ratios
cauchy_fit = cauchy_fitter(np.hstack(ratios))
cauchy_fit_85 = cauchy_fitter(np.hstack(ratios_85))
cauchy_fit_875 = cauchy_fitter(np.hstack(ratios_875))
cauchy_fit_95 = cauchy_fitter(np.hstack(ratios_95))


print(cauchy_fit_85, cauchy_fit_875, cauchy_fit, cauchy_fit_95)


# Do per channel
def cauchy_channel_fits(ratios, chunk=1):

    num_chans = vla_cube.shape[0]
    channels = np.arange(num_chans)
    chunked_channels = \
        np.array_split(channels,
                       [chunk * i for i in xrange(num_chans / chunk)])
    if chunked_channels[-1].size == 0:
        chunked_channels = chunked_channels[:-1]
    if chunked_channels[0].size == 0:
        chunked_channels = chunked_channels[1:]

    sf = np.empty(len(chunked_channels))
    sf_uncert = np.empty(len(chunked_channels))

    for i, chunk in enumerate(chunked_channels):
        slicer = slice(chunk[0], chunk[-1])

        vals = np.hstack(ratios[slicer])
        out = cauchy_fitter(vals, verbose=False)[:-1]
        sf[i] = out[0]
        sf_uncert[i] = out[1]

    return sf, sf_uncert


chunk = 10

sf_cauchy, sf_cauchy_uncert = cauchy_channel_fits(ratios, chunk=chunk)
sf_cauchy_85, sf_cauchy_uncert_85 = cauchy_channel_fits(ratios_85, chunk=chunk)
sf_cauchy_875, sf_cauchy_uncert_875 = \
    cauchy_channel_fits(ratios_875, chunk=chunk)
sf_cauchy_95, sf_cauchy_uncert_95 = cauchy_channel_fits(ratios_95, chunk=chunk)

half_chunk = chunk / 2
chunk_chans = np.arange(1, len(sf_cauchy) + 1) * chunk - half_chunk

twocolumn_figure()

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
ax[0, 0].errorbar(chunk_chans, sf_cauchy_85,
                  yerr=sf_cauchy_uncert_85, label="85 m",
                  alpha=0.5)
ax[0, 0].axhline(0, linestyle='--')
ax[0, 0].text(0, 0.4, "{}".format(np.round(beam_fwhm(85 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 0].set_ylabel(r"ln Scale Factor")
ax[0, 0].grid(True)
ax[0, 0].set_ylim([-0.3, 0.5])

ax[0, 1].errorbar(chunk_chans, sf_cauchy_875,
                  yerr=sf_cauchy_uncert_875, label="87.5 m",
                  alpha=0.5)
ax[0, 1].text(0, 0.4, "{}".format(np.round(beam_fwhm(87.5 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 1].axhline(0, linestyle='--')
ax[0, 1].grid(True)

ax[1, 0].errorbar(chunk_chans, sf_cauchy, yerr=sf_cauchy_uncert,
                  label="90 m", alpha=0.5)
ax[1, 0].axhline(0, linestyle='--')
ax[1, 0].text(0, 0.4, "{}".format(np.round(beam_fwhm(90 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 0].set_ylabel(r"ln Scale Factor")
ax[1, 0].set_xlabel("Channel")
ax[1, 0].grid(True)

ax[1, 1].errorbar(chunk_chans, sf_cauchy_95,
                  yerr=sf_cauchy_uncert_95, label="95 m",
                  alpha=0.5)
ax[1, 1].axhline(0, linestyle='--')
ax[1, 1].text(0, 0.4, "{}".format(np.round(beam_fwhm(95 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 1].set_xlabel("Channel")
ax[1, 1].grid(True)

fig.tight_layout()

fig.savefig(allfigs_path("scalefactors_cauchy_vla_gbt.png"))
fig.savefig(allfigs_path("scalefactors_cauchy_vla_gbt.pdf"))

plt.close()

# Plot comparison of fcal from methods per channel.

fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
ax[0, 0].errorbar(chunk_chans, sf_cauchy_85, label="Cauchy",
                  alpha=0.5)
ax[0, 0].errorbar(chans, np.log(sf_slope_85), label="Theil-Sen Fit",
                  alpha=0.5)
ax[0, 0].axhline(0, linestyle='--')
ax[0, 0].text(0, 0.3, "{}".format(np.round(beam_fwhm(85 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 0].set_ylabel(r"ln Scale Factor")
ax[0, 0].grid(True)
ax[0, 0].set_ylim([-1.0, 0.5])
ax[0, 0].legend(frameon=True, loc='upper right')

ax[0, 1].errorbar(chunk_chans, sf_cauchy_875, label="87.5 m",
                  alpha=0.5)
ax[0, 1].errorbar(chans, np.log(sf_slope_875), label="Sen-Theil Fit",
                  alpha=0.5)
ax[0, 1].text(0, 0.3, "{}".format(np.round(beam_fwhm(87.5 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[0, 1].axhline(0, linestyle='--')
ax[0, 1].grid(True)

ax[1, 0].errorbar(chunk_chans, sf_cauchy,
                  label="90 m", alpha=0.5)
ax[1, 0].errorbar(chans, np.log(sf_slope), label="Sen-Theil Fit",
                  alpha=0.5)
ax[1, 0].axhline(0, linestyle='--')
ax[1, 0].text(0, 0.3, "{}".format(np.round(beam_fwhm(90 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 0].set_ylabel(r"ln Scale Factor")
ax[1, 0].set_xlabel("Channel")
ax[1, 0].grid(True)

ax[1, 1].errorbar(chunk_chans, sf_cauchy_95, label="95 m",
                  alpha=0.5)
ax[1, 1].errorbar(chans, np.log(sf_slope_95), label="Sen-Theil Fit",
                  alpha=0.5)
ax[1, 1].axhline(0, linestyle='--')
ax[1, 1].text(0, 0.3, "{}".format(np.round(beam_fwhm(95 * u.m).to(u.arcmin), 1)),
              bbox={"boxstyle": "square", "facecolor": "w"})
ax[1, 1].set_xlabel("Channel")
ax[1, 1].grid(True)

fig.tight_layout()
fig.savefig(allfigs_path("scalefactors_comparison_vla_gbt.png"))
fig.savefig(allfigs_path("scalefactors_comparison_vla_gbt.pdf"))

plt.close()

onecolumn_figure()
plt.errorbar(chunk_chans, sf_cauchy_875, label="Cauchy",
             alpha=0.5)
plt.errorbar(chans, np.log(sf_slope_875), label="Sen-Theil Fit",
             alpha=0.5)
plt.axhline(0, linestyle='--')
plt.grid(True)
plt.legend(frameon=True, loc='upper left')
plt.xlabel("Channel")
plt.ylabel(r"ln Scale Factor")
plt.ylim([-1.0, 0.5])
plt.tight_layout()
plt.savefig(allfigs_path("scalefactors_comparison_vla_gbt_9.8arcmin.png"))
plt.savefig(allfigs_path("scalefactors_comparison_vla_gbt_9.8arcmin.pdf"))
plt.close()


# Finally, show the Cauchy fit for the 87.5 m beam, since this is where the
# f_cal=1 will come from.
cauchy_fit_875_limrange = cauchy_fitter(np.hstack(ratios_875[200:600]))

# Using a limited range of channels gives 0.99+/-0.24(since we are near 0)
# This is within range of just using the whole range.

out = hist(np.log(np.hstack(ratios_875)), bins='scott', normed=True)
plt.plot(np.arange(-3, 3, 0.01),
         stats.cauchy.pdf(np.arange(-3, 3, 0.01), *cauchy_fit_875[-1].params))
plt.grid(True)
plt.xlabel(r"ln I$_{\rm int}$ / I$_{\rm SD}$")
plt.tight_layout()
plt.savefig(allfigs_path("ratio_hist_vla_gbt_9.8arcmin.png"))
plt.savefig(allfigs_path("ratio_hist_vla_gbt_9.8arcmin.pdf"))
plt.close()

# The scale factor we adopt is...
print("Scale factor: {0}+/-{1}".format(np.exp(cauchy_fit_875[0]),
                                       np.abs(cauchy_fit_875[1] / cauchy_fit_875[0])))

# Open up the GBT cube and update the beam parameters
import astropy.io.fits as fits
gbt_hdu = fits.open(os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid_registered.fits"),
                    mode='update')
gbt_hdu[0].header.update(Beam(beam_fwhm(87.5 * u.m).to(u.deg)).to_header_keywords())
gbt_hdu.flush()
gbt_hdu.close()

# And the low-res version too
gbt_hdu = fits.open(os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_Tmb_14B088_spectralregrid_registered.fits"),
                    mode='update')
gbt_hdu[0].header.update(Beam(beam_fwhm(87.5 * u.m).to(u.deg)).to_header_keywords())
gbt_hdu.flush()
gbt_hdu.close()


default_figure()
