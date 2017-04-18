
'''
Define HI frequency once.
'''

import astropy.units as u

distance = 840 * u.kpc


def ang_to_phys(ang_size, distance=distance):
    '''
    Convert from angular to physical scales
    '''
    return (ang_size.to(u.rad).value * distance).to(u.pc)


hi_freq = 1.42040575177 * u.GHz

pb_lim = 0.5

co21_mass_conversion = 6.7 * (u.Msun / u.pc ** 2) / (u.K * u.km / u.s)

hi_mass_conversion = 0.019 * (u.M_sun / u.pc ** 2) / (u.K * u.km / u.s)

hi_coldens_Kkms = 1.82e18 * u.cm**-2 / (u.K * u.km / u.s)

cube_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.fits".format(pb_lim)

rotsub_cube_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.rotsub.fits".format(pb_lim)

centroidsub_cube_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.centroid_corrected.fits".format(pb_lim)

peakvelsub_cube_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.peakvels_corrected.fits".format(pb_lim)

mask_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked_source_mask" \
    ".fits".format(pb_lim)

rotsub_mask_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked_source_mask" \
    ".rotsub.fits".format(pb_lim)

centroidsub_mask_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked_source_mask" \
    ".centroid_corrected.fits".format(pb_lim)

peakvelssub_mask_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked_source_mask" \
    ".peakvels_corrected.fits".format(pb_lim)

moment0_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.mom0.fits".format(pb_lim)

moment1_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.mom1.fits".format(pb_lim)

lwidth_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.lwidth.fits".format(pb_lim)

skew_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.skewness.fits".format(pb_lim)

kurt_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.kurtosis.fits".format(pb_lim)

peakvels_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.peakvels.fits".format(pb_lim)

peaktemps_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.peaktemps.fits".format(pb_lim)

regridco_cube_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked_co21_full_regrid.fits".format(pb_lim)
