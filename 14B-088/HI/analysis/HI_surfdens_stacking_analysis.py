
'''
Analyze the outputs of HI_surfdens_stacking_feathered
'''

from pandas import DataFrame
import matplotlib.pyplot as p
import numpy as np
from spectral_cube import SpectralCube, Projection
import astropy.units as u
from astropy.io import fits

from cube_analysis.spectral_stacking_models import fit_hwhm

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, alltables_path)

from plotting_styles import default_figure, onecolumn_figure

# Compare properties of the stacked profiles
# Finally, fit Gaussian models and save the fit results

hi_mom0_hdu = fits.open(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.moment0_Kkms.fits"))[0]
hi_mom0 = Projection.from_hdu(hi_mom0_hdu)

dperc = 5
unit = hi_mom0.unit
inneredge = np.nanpercentile(hi_mom0, np.arange(0, 101, dperc)[:-1]) * unit
outeredge = np.nanpercentile(hi_mom0, np.arange(0, 101, dperc)[1:]) * unit
# Add something small to the 100th percentile so it is used
outeredge[-1] += 1e-3 * unit

wstring = "{}percentile".format(int(dperc))

sigma_noise = 2.8  # K

npix_beam = 41.

num_pix = np.load(fourteenB_HI_data_wGBT_path("stacked_spectra/surfdens_stacking_{}_num_pix.npy".format(wstring)))

peakvel_stack = SpectralCube.read(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_peak_{}.fits".format(wstring)))

hi_params = {}
hwhm_models = {}
labels = ["peaksub"]

hi_params = {}

param_names = ["sigma", "v_peak", "f_wings", "sigma_wing", "asymm", "kappa"]

for sub in labels:
    for name in param_names:
        par_name = "{0}_{1}".format(sub, name)
        par_lowlim = "{}_low_lim".format(par_name)
        par_uplim = "{}_up_lim".format(par_name)

        hi_params[par_name] = np.zeros_like(inneredge.value)
        hi_params[par_lowlim] = np.zeros_like(inneredge.value)
        hi_params[par_uplim] = np.zeros_like(inneredge.value)


for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On {0} of {1}".format(ctr + 1, len(inneredge)))

    hi_spectra = [peakvel_stack[:, ctr, 0]]

    for spectrum, label in zip(hi_spectra, labels):

        vels = spectrum.spectral_axis.to(u.km / u.s).value

        nbeams = num_pix[ctr] / npix_beam

        # Fit +/- 60 km/s
        vel_mask = np.logical_and(vels >= -100, vels <= 100)

        parvals_hwhm, parerrs_hwhm, parnames_hwhm, g_HI_hwhm, samps = \
            fit_hwhm(vels[vel_mask], spectrum.value[vel_mask],
                     sigma_noise=sigma_noise,
                     nbeams=nbeams, niters=100, interp_factor=1.)

        for idx, name in enumerate(parnames_hwhm):
            par_name = "{0}_{1}".format(label, name)
            hi_params[par_name][ctr] = parvals_hwhm[idx]
            hi_params["{}_low_lim".format(par_name)][ctr] = \
                np.abs(parerrs_hwhm[0, idx])
            hi_params["{}_up_lim".format(par_name)][ctr] = \
                np.abs(parerrs_hwhm[1, idx])

bin_names = ["{:.2f}-{:.2f} K".format(r0.value, r1.value)
             for r0, r1 in zip(inneredge, outeredge)]

bin_center = (0.5 * (inneredge + outeredge)).value
hi_params["bin_center"] = bin_center

hi_peak_fits = DataFrame(hi_params, index=bin_names)

hi_peak_fits.to_latex(alltables_path("hi_hwhm_totalprof_fits_surfdens_{}_feather.tex".format(wstring)))
hi_peak_fits.to_csv(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_surfdens_{}_feather.csv".format(wstring),
                                                no_check=True))

# Let's plot some properties.

# from pandas import read_csv
# hi_peak_fits = read_csv(fourteenB_HI_data_wGBT_path("tables/hi_hwhm_totalprof_fits_peak_{}_feather.csv".format(wstring)), index_col=0)

onecolumn_figure()

# These errorbars looked small on the plot. Unsure if it is treating it as
# single-sided or not. Doesn't really matter in this case.
p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['peaksub_sigma'],
           yerr=hi_peak_fits['peaksub_sigma_low_lim'] * 2,
           label="Peak Vel.\nSubtracted", fmt='-.^')
p.legend(frameon=True)
p.ylabel(r"$\sigma_{\rm HWHM}$ (km/s)")
p.xlabel("Integrated Intensity (K km/s)")
p.grid()
p.tight_layout()

p.savefig(allfigs_path("stacked_profiles/hi_veldisp_peak_stackedfits_feather.png"))
p.savefig(allfigs_path("stacked_profiles/hi_veldisp_peak_stackedfits_feather.pdf"))
p.close()

p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['peaksub_v_peak'],
           yerr=hi_peak_fits['peaksub_v_peak_low_lim'] * 2,
           label="Peak Vel.\nSubtracted", fmt='-.^')
p.legend(frameon=True)
p.ylabel("Centroid (km/s)")
p.xlabel("Integrated Intensity (K km/s)")
p.grid()
p.tight_layout()

p.savefig(allfigs_path("stacked_profiles/hi_vpeak_peak_stackedfits_feather.png"))
p.savefig(allfigs_path("stacked_profiles/hi_vpeak_peak_stackedfits_feather.pdf"))
p.close()

p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['peaksub_f_wings'],
           yerr=[hi_peak_fits['peaksub_f_wings_low_lim'],
                 hi_peak_fits['peaksub_f_wings_up_lim']],
           label="Peak Vel.\nSubtracted", fmt='-.^')
p.legend(frameon=True)
p.ylabel(r"$f_{\rm wings}$")
p.xlabel("Integrated Intensity (K km/s)")
p.grid()
p.tight_layout()

p.savefig(allfigs_path("stacked_profiles/hi_fwings_peak_stackedfits_feather.png"))
p.savefig(allfigs_path("stacked_profiles/hi_fwings_peak_stackedfits_feather.pdf"))
p.close()

p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['peaksub_asymm'],
           yerr=[hi_peak_fits['peaksub_asymm_low_lim'],
                 hi_peak_fits['peaksub_asymm_up_lim']],
           label="Peak Vel.\nSubtracted", fmt='-.^')
p.legend(frameon=True)
p.ylabel(r"Asymmetry")
p.xlabel("Integrated Intensity (K km/s)")
p.grid()
p.tight_layout()

p.savefig(allfigs_path("stacked_profiles/hi_asymm_peak_stackedfits_feather.png"))
p.savefig(allfigs_path("stacked_profiles/hi_asymm_peak_stackedfits_feather.pdf"))
p.close()


p.errorbar(hi_peak_fits['bin_center'], hi_peak_fits['peaksub_kappa'],
           yerr=[hi_peak_fits['peaksub_kappa_low_lim'],
                 hi_peak_fits['peaksub_kappa_up_lim']],
           label="Peak Vel.\nSubtracted", fmt='-.^')
p.legend(frameon=True)
p.ylabel(r"$\kappa$")
p.xlabel("Integrated Intensity (K km/s)")
p.grid()
p.tight_layout()

p.savefig(allfigs_path("stacked_profiles/hi_kappa_peak_stackedfits_feather.png"))
p.savefig(allfigs_path("stacked_profiles/hi_kappa_peak_stackedfits_feather.pdf"))
p.close()

default_figure()
