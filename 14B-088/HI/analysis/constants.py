
'''
Define HI frequency once.
'''

import astropy.units as u

hi_freq = 1.42040575177 * u.GHz

pb_lim = 0.5

cube_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.fits".format(pb_lim)

rotsub_cube_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked.rotsub.fits".format(pb_lim)

mask_name = \
    "M33_14B-088_HI.clean.image.pbcov_gt_{}_masked_source_mask.fits".format(pb_lim)

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
