from spectral_cube import SpectralCube, Projection
import numpy as np
import matplotlib.pyplot as p
from astropy.io import fits
import astropy.units as u
from corner import hist2d

from paths import (fourteenB_HI_data_path, fourteenB_HI_file_dict,
                   fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, iram_co21_data_path)
from constants import hi_freq
from plotting_styles import (twocolumn_twopanel_figure, onecolumn_figure,
                             default_figure)

'''
Investigating skewness and kurtosis in the 14B-088 cube.
'''

cube = SpectralCube.read(fourteenB_HI_file_dict["Cube"])
mask = fits.open(fourteenB_HI_file_dict['Source_Mask'])[0].data > 0
cube = cube.with_mask(mask)

mom0_hdu = fits.open(fourteenB_HI_file_dict["Moment0"])[0]
mom0 = Projection.from_hdu(mom0_hdu)

mom1_hdu = fits.open(fourteenB_HI_file_dict["Moment1"])[0]
mom1 = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_HI_file_dict["LWidth"])[0]
lwidth = Projection.from_hdu(lwidth_hdu)

skew_hdu = fits.open(fourteenB_HI_file_dict["Skewness"])[0]
skew = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_HI_file_dict["Kurtosis"])[0]
kurt = Projection.from_hdu(kurt_hdu)

peaktemps_hdu = fits.open(fourteenB_HI_file_dict["PeakTemp"])[0]
peaktemps = Projection.from_hdu(peaktemps_hdu)

peakvels_hdu = fits.open(fourteenB_HI_file_dict["PeakVels"])[0]
peakvels = Projection.from_hdu(peakvels_hdu)


# Compare the population of values to different properties
# Convert zeroth moment to K km/s from Jy m/s
mom0_vals = (mom0.array.flatten() * mom0.beam.jtok(hi_freq) *
             u.km / u.s / 1000.).value
mom1_vals = mom1.to(u.km / u.s).array.flatten()
lwidth_vals = lwidth.to(u.km / u.s).array.flatten()
skew_vals = skew.array.flatten()
kurt_vals = kurt.array.flatten()
peaktemps_vals = peaktemps.array.flatten()
peakvels_vals = peakvels.to(u.km / u.s).array.flatten()

good_vals = np.isfinite(skew_vals)

mask_summed = cube.mask.include().sum(0)

# Feathered versions
cube_feath = SpectralCube.read(fourteenB_wGBT_HI_file_dict["Cube"])
mask_feath = fits.open(fourteenB_wGBT_HI_file_dict['Source_Mask'])[0].data > 0
cube_feath = cube_feath.with_mask(mask_feath)

mom0_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0]
mom0_feath = Projection.from_hdu(mom0_hdu)

mom1_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Moment1"])[0]
mom1_feath = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_wGBT_HI_file_dict["LWidth"])[0]
lwidth_feath = Projection.from_hdu(lwidth_hdu)

skew_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Skewness"])[0]
skew_feath = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_wGBT_HI_file_dict["Kurtosis"])[0]
kurt_feath = Projection.from_hdu(kurt_hdu)

peaktemps_hdu = fits.open(fourteenB_wGBT_HI_file_dict["PeakTemp"])[0]
peaktemps_feath = Projection.from_hdu(peaktemps_hdu)

peakvels_hdu = fits.open(fourteenB_wGBT_HI_file_dict["PeakVels"])[0]
peakvels_feath = Projection.from_hdu(peakvels_hdu)


# Compare the population of values to different properties
# Convert zeroth moment to K km/s from Jy m/s
mom0_feath_vals = (mom0_feath.array.flatten() * mom0.beam.jtok(hi_freq) *
                   u.km / u.s / 1000.).value
mom1_feath_vals = mom1_feath.to(u.km / u.s).array.flatten()
lwidth_feath_vals = lwidth_feath.to(u.km / u.s).array.flatten()
skew_feath_vals = skew_feath.array.flatten()
kurt_feath_vals = kurt_feath.array.flatten()
peaktemps_feath_vals = peaktemps_feath.array.flatten()
peakvels_feath_vals = peakvels_feath.to(u.km / u.s).array.flatten()

good_feath_vals = np.isfinite(skew_feath_vals)

mask_summed_feath = cube_feath.mask.include().sum(0)

# onecolumn_figure(font_scale=1.1)
twocolumn_twopanel_figure(font_scale=1.1)
# Skew vs. kurt
fig, ax = p.subplots(1, 2)
hist2d(skew_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-3, 3), (-3, 6)], ax=ax[0])
ax[0].set_xlabel("Skewness")
ax[0].set_ylabel("Kurtosis")
hist2d(skew_feath_vals[good_feath_vals], kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-3, 3), (-3, 6)], ax=ax[1])
ax[1].set_xlabel("Skewness")
# ax[1].set_ylabel("Kurtosis")
p.tight_layout()

p.savefig(allfigs_path("skewness_vs_kurtosis.pdf"))
p.savefig(allfigs_path("skewness_vs_kurtosis.png"))
p.close()

# Versus integrated intensity
fig, ax = p.subplots(1, 2)
hist2d(mom0_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_vals)), (-3, 3)],
       ax=ax[0])
ax[0].set_xlabel(r"Integrated Intensity (K km s$^{-1}$)")
ax[0].set_ylabel("Skewness")
hist2d(mom0_feath_vals[good_feath_vals], skew_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_feath_vals)), (-3, 3)],
       ax=ax[1])
ax[1].set_xlabel(r"Integrated Intensity (K km s$^{-1}$)")
p.tight_layout()

p.savefig(allfigs_path("mom0_vs_skewness.pdf"))
p.savefig(allfigs_path("mom0_vs_skewness.png"))
p.close()

fig, ax = p.subplots(1, 2)
hist2d(mom0_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_vals)), (-3, 6)],
       ax=ax[0])
ax[0].set_xlabel(r"Integrated Intensity (K km s$^{-1}$)")
ax[0].set_ylabel("Kurtosis")
hist2d(mom0_feath_vals[good_feath_vals], kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_feath_vals)), (-3, 6)],
       ax=ax[1])
ax[1].set_xlabel(r"Integrated Intensity (K km s$^{-1}$)")
p.tight_layout()

p.savefig(allfigs_path("mom0_vs_kurtosis.pdf"))
p.savefig(allfigs_path("mom0_vs_kurtosis.png"))
p.close()

# Take a bin of summed mask values and see where this lies on the mom0 vs
# kurtosis plane
fig, ax = p.subplots(1, 2)
hist2d(mom0_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_vals)), (-3, 6)],
       ax=ax[0])
ax[0].set_xlabel(r"Integrated Intensity (K km s$^{-1}$)")
ax[0].set_ylabel("Kurtosis")
hist2d(mom0_feath_vals[good_feath_vals], kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(-1, np.nanmax(mom0_feath_vals)), (-3, 6)],
       ax=ax[1])
ax[1].set_xlabel(r"Integrated Intensity (K km s$^{-1}$)")
p.tight_layout()

mask_summed_300 = mask_summed.flatten()[good_vals] == 300
ax[0].plot(mom0_vals[good_vals][mask_summed_300],
           kurt_vals[good_vals][mask_summed_300], 'D')
mask_summed_300_feath = mask_summed_feath.flatten()[good_feath_vals] == 300
ax[1].plot(mom0_feath_vals[good_feath_vals][mask_summed_300_feath],
           kurt_feath_vals[good_feath_vals][mask_summed_300_feath], 'D')
p.tight_layout()
p.savefig(allfigs_path("mom0_vs_kurtosis_w_maskwidth300.pdf"))
p.savefig(allfigs_path("mom0_vs_kurtosis_w_maskwidth300.png"))
p.close()

# Compare velocity surfaces to skew and kurt.
fig, ax = p.subplots(1, 2)
hist2d(mom1_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mom1_vals), np.nanmax(mom1_vals)), (-3, 3)],
       ax=ax[0])
ax[0].set_xlabel(r"Centroid Velocity (km s$^{-1}$)")
ax[0].set_ylabel("Skewness")
hist2d(mom1_feath_vals[good_feath_vals], skew_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mom1_feath_vals),
               np.nanmax(mom1_feath_vals)), (-3, 3)],
       ax=ax[1])
ax[1].set_xlabel(r"Centroid Velocity (km s$^{-1}$)")
p.tight_layout()

p.savefig(allfigs_path("mom1_vs_skewness.pdf"))
p.savefig(allfigs_path("mom1_vs_skewness.png"))
p.close()

fig, ax = p.subplots(1, 2)
hist2d(mom1_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mom1_vals), np.nanmax(mom1_vals)), (-5, 9)],
       ax=ax[0])
ax[0].set_xlabel(r"Centroid Velocity (km s$^{-1}$)")
ax[0].set_ylabel("Kurtosis")
hist2d(mom1_feath_vals[good_feath_vals], kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mom1_feath_vals),
               np.nanmax(mom1_feath_vals)), (-5, 9)],
       ax=ax[1])
ax[1].set_xlabel(r"Centroid Velocity (km s$^{-1}$)")
p.tight_layout()

p.savefig(allfigs_path("mom1_vs_kurtosis.pdf"))
p.savefig(allfigs_path("mom1_vs_kurtosis.png"))
p.close()

fig, ax = p.subplots(1, 2)
hist2d(peakvels_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peakvels_vals), np.nanmax(peakvels_vals)), (-3, 3)],
       ax=ax[0])
ax[0].set_xlabel(r"Peak Velocity (km s$^{-1}$)")
ax[0].set_ylabel("Skewness")
hist2d(peakvels_feath_vals[good_feath_vals], skew_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peakvels_feath_vals),
               np.nanmax(peakvels_feath_vals)), (-3, 3)],
       ax=ax[1])
ax[1].set_xlabel(r"Peak Velocity (km s$^{-1}$)")
p.tight_layout()

p.savefig(allfigs_path("peakvel_vs_skewness.pdf"))
p.savefig(allfigs_path("peakvel_vs_skewness.png"))
p.close()


fig, ax = p.subplots(1, 2)
hist2d(peakvels_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peakvels_vals), np.nanmax(peakvels_vals)), (-5, 9)],
       ax=ax[0])
ax[0].set_xlabel(r"Peak Velocity (km s$^{-1}$)")
ax[0].set_ylabel("Kurtosis")
hist2d(peakvels_feath_vals[good_feath_vals], kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peakvels_feath_vals),
               np.nanmax(peakvels_feath_vals)), (-5, 9)],
       ax=ax[1])
ax[1].set_xlabel(r"Peak Velocity (km s$^{-1}$)")
p.tight_layout()

p.savefig(allfigs_path("peakvel_vs_kurtosis.pdf"))
p.savefig(allfigs_path("peakvel_vs_kurtosis.png"))
p.close()

# Versus peak temperature
fig, ax = p.subplots(1, 2)
hist2d(peaktemps_vals[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps_vals), np.nanmax(peaktemps_vals)), (-3, 3)],
       ax=ax[0])
ax[0].set_xlabel("Peak Temperature (K)")
ax[0].set_ylabel("Skewness")
hist2d(peaktemps_feath_vals[good_feath_vals], skew_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps_feath_vals),
               np.nanmax(peaktemps_feath_vals)), (-3, 3)],
       ax=ax[1])
ax[1].set_xlabel("Peak Temperature (K)")
p.tight_layout()

p.savefig(allfigs_path("peaktemp_vs_skewness.pdf"))
p.savefig(allfigs_path("peaktemp_vs_skewness.png"))
p.close()

fig, ax = p.subplots(1, 2)
hist2d(peaktemps_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps_vals), np.nanmax(peaktemps_vals)), (-5, 9)],
       ax=ax[0])
ax[0].set_xlabel("Peak Temperature (K)")
ax[0].set_ylabel("Kurtosis")
hist2d(peaktemps_feath_vals[good_feath_vals], kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps_feath_vals),
               np.nanmax(peaktemps_feath_vals)), (-5, 9)],
       ax=ax[1])
ax[1].set_xlabel("Peak Temperature (K)")
p.tight_layout()

p.savefig(allfigs_path("peaktemp_vs_kurtosis.pdf"))
p.savefig(allfigs_path("peaktemp_vs_kurtosis.png"))
p.close()

fig, ax = p.subplots(1, 2)
hist2d(peaktemps_vals[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps_vals), np.nanmax(peaktemps_vals)), (-5, 9)],
       ax=ax[0])
ax[0].set_xlabel("Peak Temperature (K)")
ax[0].set_ylabel("Kurtosis")
hist2d(peaktemps_feath_vals[good_feath_vals], kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(peaktemps_feath_vals),
               np.nanmax(peaktemps_feath_vals)), (-5, 9)],
       ax=ax[1])
ax[1].set_xlabel("Peak Temperature (K)")
p.tight_layout()

ax[0].plot(peaktemps_vals[good_vals][mask_summed_300],
           kurt_vals[good_vals][mask_summed_300], 'D')
ax[1].plot(peaktemps_feath_vals[good_feath_vals][mask_summed_300_feath],
           kurt_feath_vals[good_feath_vals][mask_summed_300_feath], 'D')

p.tight_layout()
p.savefig(allfigs_path("peaktemp_vs_kurtosis_w_maskwidth300.pdf"))
p.savefig(allfigs_path("peaktemp_vs_kurtosis_w_maskwidth300.png"))
p.close()


# How badly is each affected by the shape of the mask?
fig, ax = p.subplots(1, 2)
hist2d(mask_summed.flatten()[good_vals], skew_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)), (-3, 3)],
       ax=ax[0])
ax[0].set_xlabel("Summed Mask (spectral pixels)")
ax[0].set_ylabel("Skewness")
hist2d(mask_summed_feath.flatten()[good_feath_vals],
       skew_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mask_summed_feath), np.nanmax(mask_summed_feath)),
              (-3, 3)],
       ax=ax[1])
ax[1].set_xlabel("Summed Mask (spectral pixels)")
p.tight_layout()

p.savefig(allfigs_path("maskshape_vs_skewness.pdf"))
p.savefig(allfigs_path("maskshape_vs_skewness.png"))
p.close()

fig, ax = p.subplots(1, 2)
hist2d(mask_summed.flatten()[good_vals], kurt_vals[good_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)), (-5, 9)],
       ax=ax[0])
ax[0].set_xlabel("Summed Mask (spectral pixels)")
ax[0].set_ylabel("Kurtosis")
hist2d(mask_summed_feath.flatten()[good_feath_vals],
       kurt_feath_vals[good_feath_vals],
       data_kwargs={"alpha": 0.2},
       range=[(np.nanmin(mask_summed_feath), np.nanmax(mask_summed_feath)),
              (-5, 9)],
       ax=ax[1])
ax[1].set_xlabel("Summed Mask (spectral pixels)")
p.tight_layout()


p.savefig(allfigs_path("maskshape_vs_kurtosis.pdf"))
p.savefig(allfigs_path("maskshape_vs_kurtosis.png"))
p.close()

# fig, ax = p.subplots(1, 2)
# hist2d(mask_summed.flatten()[good_vals], peaktemps_vals[good_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)),
#               (np.nanmin(peaktemps_vals), np.nanmax(peaktemps_vals))],
#        ax=ax[0])
# ax[0].set_xlabel("Summed Mask (spectral pixels)")
# ax[0].set_ylabel("Peak Temperature")
# hist2d(mask_summed_feath.flatten()[good_feath_vals],
#        peaktemps_feath_vals[good_feath_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mask_summed_feath), np.nanmax(mask_summed_feath)),
#               (np.nanmin(peaktemps_feath_vals),
#                np.nanmax(peaktemps_feath_vals))],
#        ax=ax[1])
# ax[1].set_xlabel("Summed Mask (spectral pixels)")

# fig, ax = p.subplots(1, 2, sharey=True, sharex=True)
# hist2d(mask_summed.flatten()[good_vals], mom0_vals[good_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mask_summed), np.nanmax(mask_summed)),
#               (np.nanmin(mom0_vals), np.nanmax(mom0_vals))],
#        ax=ax[0])
# ax[0].set_xlabel("Summed Mask (spectral pixels)")
# ax[0].set_ylabel(r"Integrated Intensity (K km s$^{-1}$)")
# hist2d(mask_summed_feath.flatten()[good_feath_vals],
#        mom0_feath_vals[good_feath_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mask_summed_feath), np.nanmax(mask_summed_feath)),
#               (np.nanmin(mom0_feath_vals),
#                np.nanmax(mom0_feath_vals))],
#        ax=ax[1])
# ax[1].set_xlabel("Summed Mask (spectral pixels)")

# hist2d(mom0_vals[good_vals], peaktemps_vals[good_vals],
#        data_kwargs={"alpha": 0.2},
#        range=[(np.nanmin(mom0_vals), np.nanmax(mom0_vals)),
#               (np.nanmin(peaktemps_vals), np.nanmax(peaktemps_vals))])

default_figure()
