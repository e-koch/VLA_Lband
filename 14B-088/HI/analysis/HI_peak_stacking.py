
from spectral_cube import SpectralCube, Projection
import numpy as np
import astropy.units as u
from astropy.io import fits

from astropy.modeling import models, fitting
from pandas import DataFrame
import matplotlib.pyplot as p
import os

from cube_analysis.spectral_stacking import total_profile

from paths import (fourteenB_HI_data_path, fourteenB_HI_file_dict,
                   allfigs_path, alltables_path)

from constants import hi_freq

from plotting_styles import default_figure, onecolumn_figure


hi_cube = SpectralCube.read(fourteenB_HI_file_dict["RotSub_Cube"])
hi_mask = fits.open(fourteenB_HI_file_dict["RotSub_Mask"])[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

hi_cube_cent = SpectralCube.read(fourteenB_HI_file_dict["CentSub_Cube"])
hi_mask_cent = fits.open(fourteenB_HI_file_dict["CentSub_Mask"])[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

hi_cube_peakvel = SpectralCube.read(fourteenB_HI_file_dict["PeakSub_Cube"])
hi_mask_peakvel = fits.open(fourteenB_HI_file_dict["PeakSub_Mask"])[0]
hi_cube_peakvel = hi_cube_cent.with_mask(hi_mask_peakvel.data > 0)

hi_peaktemp_hdu = fits.open(fourteenB_HI_file_dict["PeakTemp"])[0]
hi_peaktemp = Projection.from_hdu(hi_peaktemp_hdu)

hi_beam = hi_cube.beam

dperc = 5
unit = hi_peaktemp.unit
inneredge = np.nanpercentile(hi_peaktemp, np.arange(0, 101, dperc)[:-1]) * unit
outeredge = np.nanpercentile(hi_peaktemp, np.arange(0, 101, dperc)[1:]) * unit
# Add something small to the 100th percentile so it is used
outeredge[-1] += 1e-3 * unit

total_spectrum_hi_peak = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_hi_peak_cent = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_hi_peak_peakvel = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K

for ctr, (p0, p1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {} K".format(p0, p1))

    mask = np.logical_and(hi_peaktemp >= p0, hi_peaktemp < p1)

    total_spectrum_hi_peak[ctr] = \
        total_profile(hi_cube, mask, num_cores=4).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_peak_cent[ctr] = \
        total_profile(hi_cube_cent, mask, num_cores=4).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_peak_peakvel[ctr] = \
        total_profile(hi_cube_peakvel, mask, num_cores=4).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))


# We'll make mock SpectralCubes from these so it's easy to calculate
# moments and such from
rot_stack = SpectralCube(data=total_spectrum_hi_peak.T.reshape((1178, inneredge.size, 1)),
                         wcs=hi_cube.wcs)

cent_stack = SpectralCube(data=total_spectrum_hi_peak_cent.T.reshape((1178, inneredge.size, 1)),
                          wcs=hi_cube.wcs)

peakvel_stack = SpectralCube(data=total_spectrum_hi_peak_peakvel.T.reshape((1178, inneredge.size, 1)),
                             wcs=hi_cube.wcs)

# Now save all of these for future use.
stacked_folder = fourteenB_HI_data_path("stacked_spectra", no_check=True)
if os.path.exists(stacked_folder):
    os.mkdir(stacked_folder)

wstring = "{}percentile".format(int(dperc))
rot_stack.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_peak_{}.fits".format(wstring),
                                       no_check=True), overwrite=True)

cent_stack.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_peak_{}.fits".format(wstring),
                                        no_check=True), overwrite=True)

peakvel_stack.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_peak_{}.fits".format(wstring),
                                           no_check=True), overwrite=True)


# Compare properties of the stacked profiles
# Finally, fit Gaussian models and save the fit results

g_HI_init = models.Gaussian1D(amplitude=1., mean=0., stddev=10.)

hi_params = {}
labels = ["rotsub", "centsub", "peaksub"]

for sub in labels:
    for name in g_HI_init.param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_error = "{}_stderr".format(par_name)

        hi_params[par_name] = np.zeros_like(inneredge.value)
        hi_params[par_error] = np.zeros_like(inneredge.value)

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    hi_spectra = [rot_stack[:, ctr, 0],
                  cent_stack[:, ctr, 0],
                  peakvel_stack[:, ctr, 0]]

    for spectrum, label in zip(hi_spectra, labels):

        fit_g = fitting.LevMarLSQFitter()

        vels = hi_cube.spectral_axis.to(u.km / u.s).value
        norm_intens = (spectrum / spectrum.max()).value
        g_HI = fit_g(g_HI_init, vels, norm_intens, maxiter=1000)

        cov = fit_g.fit_info['param_cov']
        if cov is None:
            raise Exception("No covariance matrix")

        idx_corr = 0
        for idx, name in enumerate(g_HI.param_names):
            if name == "mean_1":
                idx_corr = 1
                continue
            par_name = "{0}_{1}".format(label, name)
            hi_params[par_name][ctr] = g_HI.parameters[idx]
            hi_params["{}_stderr".format(par_name)][ctr] = \
                np.sqrt(cov[idx - idx_corr, idx - idx_corr])

bin_names = ["{:.2f}-{:.2f} K".format(r0.value, r1.value)
             for r0, r1 in zip(inneredge, outeredge)]

bin_center = (0.5 * (inneredge + outeredge)).value
hi_params["bin_center"] = bin_center

# Add stderr in quadrature with the channel width
hi_velres = \
    (hi_cube.spectral_axis[1] -
     hi_cube.spectral_axis[0]).to(u.km / u.s).value

# Add the velocity width of the channel in quadrature
for col in hi_params.keys():
    if "mean_stderr" in col or "stddev_stderr" in col:
        hi_params[col + "_w_chanwidth"] = np.sqrt(hi_params[col]**2 +
                                                  hi_velres**2)

hi_peak_fits = DataFrame(hi_params, index=bin_names)

hi_peak_fits.to_latex(alltables_path("hi_gaussian_totalprof_fits_peak_{}.tex".format(wstring)))
hi_peak_fits.to_csv(fourteenB_HI_data_path("tables/hi_gaussian_totalprof_fits_peak_{}.csv".format(wstring),
                                           no_check=True))

# Let's plot some properties.

onecolumn_figure()

p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['rotsub_stddev'],
           yerr=hi_peak_fits['rotsub_stddev_stderr_w_chanwidth'],
           label="Rotation\nSubtracted")
p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['centsub_stddev'],
           yerr=hi_peak_fits['centsub_stddev_stderr_w_chanwidth'],
           label="Centroid\nSubtracted")
p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['peaksub_stddev'],
           yerr=hi_peak_fits['peaksub_stddev_stderr_w_chanwidth'],
           label="Peak Vel.\nSubtracted")
p.legend(frameon=True)
p.ylabel("Velocity Dispersion (km/s)")
p.xlabel("Peak Temperature (K)")
p.grid()
p.tight_layout()

p.savefig(allfigs_path("hi_veldisp_peak_stackedfits.png"))
p.savefig(allfigs_path("hi_veldisp_peak_stackedfits.pdf"))
p.close()

p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['rotsub_mean'],
           yerr=hi_peak_fits['rotsub_mean_stderr_w_chanwidth'],
           label="Rotation\nSubtracted")
p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['centsub_mean'],
           yerr=hi_peak_fits['centsub_mean_stderr_w_chanwidth'],
           label="Centroid\nSubtracted")
p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['peaksub_mean'],
           yerr=hi_peak_fits['peaksub_mean_stderr_w_chanwidth'],
           label="Peak Vel.\nSubtracted")
p.legend(frameon=True)
p.ylabel("Centroid (km/s)")
p.xlabel("Peak Temperature (K)")
p.grid()
p.tight_layout()

p.savefig(allfigs_path("hi_centroid_peak_stackedfits.png"))
p.savefig(allfigs_path("hi_centroid_peak_stackedfits.pdf"))
p.close()

default_figure()
