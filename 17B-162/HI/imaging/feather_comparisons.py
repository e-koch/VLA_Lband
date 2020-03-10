
'''
Compare the data where they overlap in the uv plane.

No offset correction is needed.
'''

from spectral_cube import SpectralCube
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import os
import scipy.ndimage as nd

from uvcombine.scale_factor import find_scale_factor

from cube_analysis.feather_cubes import feather_compare_cube

from paths import (seventeenB_HI_data_02kms_path,
                   seventeenB_HI_data_1kms_path,
                   data_path, allfigs_path)
from constants import hi_freq
from plotting_styles import onecolumn_figure

# Compare with the 1 km/s cube. Higher S/N
# vla_cube = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.fits"))
vla_cube = SpectralCube.read(seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.image.pbcor.fits"))

# pb_cube = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.pb.fits"))
pb_cube = SpectralCube.read(seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.pb.fits"))
# PB minimally changes over the frequency range. So just grab one plane
pb_plane = pb_cube[0]

# We need to define a tapered weighting function to ignore emission outside
# of the VLA mosaic


def taper_weights(mask, sigma, nsig_cut=3):

    dist = nd.distance_transform_edt(mask)

    gauss_dists = np.where(np.logical_and(dist < nsig_cut * sigma, dist > 0.))
    flat_dists = np.where(dist >= nsig_cut * sigma)

    weight_arr = np.zeros_like(mask, dtype=float)

    weight_arr[gauss_dists] = \
        np.exp(- (dist[gauss_dists] - nsig_cut * sigma)**2 / (2 * sigma**2))
    weight_arr[flat_dists] = 1.

    return weight_arr


weight = taper_weights(np.isfinite(pb_plane), 30, nsig_cut=5)


gbt_path = os.path.join(data_path, "GBT")
# gbt_cube = SpectralCube.read(os.path.join(gbt_path, "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_02kms.fits"))
gbt_cube = SpectralCube.read(os.path.join(gbt_path, "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms.fits"))

beam_fwhm = lambda diam: ((1.18 * hi_freq.to(u.cm, u.spectral())) / diam.to(u.cm)) * u.rad

# Already determined from the 14B HI analysis. Lowered spatial resolution
# due to lack of overlap in the GBT fields centered at M33. So the data were
# gridded with a Gaussian kernel, rather than a jinc function
gbt_eff_beam = beam_fwhm(87.5 * u.m)

# The shortest baseline in the 14B-088 data is ~44 m.
las = (hi_freq.to(u.cm, u.spectral()) / (44 * u.m)).to(u.arcsec, u.dimensionless_angles())

radii, ratios, high_pts, low_pts, chan_out = \
    feather_compare_cube(vla_cube, gbt_cube, las,
                         num_cores=1,
                         lowresfwhm=gbt_eff_beam,
                         chunk=50,
                         verbose=False,
                         weights=weight,
                         relax_spectral_check=False,
                         # NOTE: there is an offset of ~0.4 km/s between the cubes
                         # The big GBT beam means this really doesn't matter (I
                         # manually checked). The difference is 0.36 times the
                         # channel size. I have no idea where this shift is coming
                         # from since the freq axis used in `gbt_regrid.py` matches
                         # the frequency in the individual channel MSs used in
                         # imaging. It's not even a half-channel offset like I
                         # would expect if the MS frequency was the channel edge...
                         spec_check_kwargs={'rtol': 0.4})

onecolumn_figure()
sc_factor, sc_err = find_scale_factor(np.hstack(low_pts), np.hstack(high_pts),
                                      method='distrib',
                                      verbose=True)

plt.grid(True)
plt.xlabel(r"ln I$_{\rm int}$ / I$_{\rm SD}$")
plt.tight_layout()
plt.savefig(allfigs_path("Imaging/ratio_hist_17B_vla_gbt_9.8arcmin_v3_w_weights.png"))
plt.savefig(allfigs_path("Imaging/ratio_hist_17B_vla_gbt_9.8arcmin_v3_w_weights.pdf"))

print("Factor: {0}+/-{1}".format(sc_factor, sc_err))
# Factor: 1.125046+/-0.00394768

# This isn't a fantastic fit, so this error was significantly underestimated

plt.close()

# Compare properties per-channel
sc_factor_chans = []
sc_err_chans = []
for low, high in zip(low_pts, high_pts):
    sc_f, sc_e = \
        find_scale_factor(low, high,
                          method='distrib',
                          verbose=False)
    sc_factor_chans.append(sc_f)
    sc_err_chans.append(sc_e)


sc_factor_chans_linfit = []
sc_err_chans_linfit = []
for low, high in zip(low_pts, high_pts):
    sc_f, sc_e = \
        find_scale_factor(low, high,
                          method='linfit',
                          verbose=False)
    sc_factor_chans_linfit.append(sc_f)
    sc_err_chans_linfit.append(sc_e)

sc_factor_chans_linfit = np.array(sc_factor_chans_linfit)
sc_err_chans_linfit = np.array(sc_err_chans_linfit)

chans = np.arange(len(low_pts))

onecolumn_figure()
plt.errorbar(chans, sc_factor_chans,
             yerr=sc_err_chans,
             alpha=0.5, label='Distrib Fit')
plt.errorbar(chans, sc_factor_chans_linfit,
             yerr=[sc_factor_chans_linfit - sc_err_chans_linfit[:, 0],
                   sc_err_chans_linfit[:, 1] - sc_factor_chans_linfit],
             alpha=0.5, label='Linear fit')
# plt.plot(chans, slope_lowess_85)
plt.axhline(1, linestyle='--')
plt.legend(frameon=True)
plt.ylabel(r"Scale Factor")
plt.xlabel("Channels")
plt.grid(True)

plt.tight_layout()

plt.savefig(allfigs_path("Imaging/ratio_hist_perchan_17B_vla_gbt_9.8arcmin_v3_w_weights.png"))
plt.savefig(allfigs_path("Imaging/ratio_hist_perchan_17B_vla_gbt_9.8arcmin_v3_w_weights.pdf"))
plt.close()

# Now refit with the channels near the systemic velocity, where most of the HI
# structure falls within the mosaic PB
chan_range = slice(80, 160)

onecolumn_figure()
sc_factor_chrange, sc_err_chrange = \
    find_scale_factor(np.hstack(low_pts[chan_range]),
                      np.hstack(high_pts[chan_range]),
                      method='distrib',
                      verbose=True)

plt.grid(True)
plt.xlabel(r"ln I$_{\rm int}$ / I$_{\rm SD}$")
plt.tight_layout()
plt.savefig(allfigs_path(f"Imaging/ratio_hist_17B_vla_gbt_9.8arcmin_chan_{chan_range.start}_{chan_range.stop}_v3_w_weights.png"))
plt.savefig(allfigs_path(f"Imaging/ratio_hist_17B_vla_gbt_9.8arcmin_chan_{chan_range.start}_{chan_range.stop}_v3_w_weights.pdf"))

print("Factor: {0}+/-{1}".format(sc_factor_chrange, sc_err_chrange))
# Factor: 1.105133+/-0.00463
# Error still underestimated

# The >1 factor is due to some emission in the GBT data being cut-off by the
# PB limit of the VLA mosaic. The factor increases far from the systemic
# velocity, where bright HI gets cut-off (compared to the larger 14B data).
# So, despite the != 1 factor, no factor will be applied to the SD data.
# Besides, the 14B mosaic comparison gives a 1.0 factor with the GBT data.
# The tests here were for consistency and that's what we find.

plt.close()
