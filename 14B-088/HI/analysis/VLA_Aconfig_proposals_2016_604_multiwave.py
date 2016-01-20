
'''
Overlaying VLA, Hubble, Chandra data

Views to examine:
 * HI, Halpha, 2-10 keV
 * HI, [SII], 2-10 keV
 * L-band continuum
 * CO(2-1)

'''

from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.wcs_utils import drop_axis

from astropy.io import fits
import astropy.units as u
from wcsaxes import WCS
from astropy.visualization import scale_image

import sys
import numpy as np
import matplotlib.pyplot as p

# Set which plot to create
# whichplot = 4
whichplot = int(sys.argv[1])


def load_hi():
    hi_mom0 = \
        fits.open("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"
                  "M33_14B-088_HI.clean.image.mom0_3sigmask.fits")[0]

    hi_proj = Projection(hi_mom0.data, unit=u.Jy, wcs=WCS(hi_mom0.header))

    return hi_mom0, hi_proj


def load_halp_hubble():

    halpha = fits.open("/media/eric/Data_3/M33/Hubble/hst_05773_03_wfpc2_f656n_wf"
                       "/hst_05773_03_wfpc2_f656n_wf_sci.fits")[1]

    halpha_proj = Projection(halpha.data, unit=u.count,
                             wcs=WCS(halpha.header))

    return halpha, halpha_proj


def load_sii_hubble():

    sii = fits.open("/media/eric/Data_3/M33/Hubble/hst_05773_03_wfpc2_f673n_wf/"
                    "hst_05773_03_wfpc2_f673n_wf_sci.fits")[1]

    sii_proj = Projection(sii.data, unit=u.count, wcs=WCS(sii.header))

    return sii, sii_proj


def load_sii_halp_hubble():

    sii_halp = fits.open("/media/eric/Data_3/M33/Hubble/"
                    "hst_05773_wfpc2_sii_halp_ratio.fits")[0]

    sii_halp_proj = Projection(sii_halp.data, unit=u.count,
                               wcs=WCS(sii_halp.header))

    return sii_halp, sii_halp_proj


def load_chandra():
    chandra = fits.open("/media/eric/Data_3/M33/Chandra/images/"
                        "broad_thresh_exp_sm.img")[0]

    chandra_proj = Projection(chandra.data, unit=u.count,
                              wcs=WCS(chandra.header))

    return chandra, chandra_proj


def load_co21_IRAM():

    co21 = fits.open("/media/eric/Data_3/M33/co21/m33.ico.fits")[0]

    new_wcs = drop_axis(drop_axis(WCS(co21.header), 3), 2)
    co21_proj = Projection(co21.data.squeeze(), unit=u.K, wcs=new_wcs)

    return co21, co21_proj


def load_halp_kp():

    halp_kp = fits.open("/media/eric/Data_3/M33/Halpha/ha.fits")[0]

    halp_kp_proj = Projection(halp_kp.data, unit=u.count,
                              wcs=WCS(halp_kp.header))

    return halp_kp, halp_kp_proj


def load_continuum_L():

    continuum = fits.open("/media/eric/Data_2/14B-088_continuum/test_imaging/"
                          "NGC_604_M33_3_briggs_0.5.image.tt0.fits")[0]

    new_wcs = drop_axis(drop_axis(WCS(continuum.header), 3), 2)

    continuum_proj = Projection(continuum.data.squeeze(), unit=u.Jy,
                                wcs=new_wcs)

    return continuum, continuum_proj


# HI, Halpha, 2-10 keV
if whichplot == 0:

    hi_mom0, hi_proj = load_hi()
    halpha, halpha_proj = load_halp_hubble()
    chandra, chandra_proj = load_chandra()

    halpha_proj.quicklook()
    halpha_proj.FITSFigure.show_contour(hi_mom0, colors='w', alpha=0.8)
    halpha_proj.FITSFigure.show_contour(chandra, colors='r', alpha=0.8)

# HI, [SII], 2-10 keV
elif whichplot == 1:

    hi_mom0, hi_proj = load_hi()
    sii, sii_proj = load_sii_hubble()
    chandra, chandra_proj = load_chandra()

    sii_proj.quicklook()
    sii_proj.FITSFigure.show_contour(hi_mom0, colors='w', alpha=0.8)
    sii_proj.FITSFigure.show_contour(chandra, colors='r', alpha=0.8, )

# [SII]/Halp, 2-10 keV, CO21
elif whichplot == 2:

    sii_halp, sii_halp_proj = load_sii_halp_hubble()
    chandra, chandra_proj = load_chandra()
    co21, co21_proj = load_co21_IRAM()

    sii_halp_proj.quicklook()
    sii_halp_proj.FITSFigure.show_contour(co21, colors='r', alpha=0.8)
    sii_halp_proj.FITSFigure.show_contour(chandra, colors='b', alpha=0.8, )

# HI, Halpha SD, CO21
elif whichplot == 3:

    hi_mom0, hi_proj = load_hi()
    halp_kp, halp_kp_proj = load_halp_kp()
    co21, co21_proj = load_co21_IRAM()

    halp_kp_proj.quicklook()
    halp_kp_proj.FITSFigure.show_contour(hi_mom0, colors='g', alpha=0.8)
    halp_kp_proj.FITSFigure.show_contour(co21, colors='r', alpha=0.8)

# HI, Halpha, CO21
elif whichplot == 4:

    halpha, halpha_proj = load_halp_hubble()
    hi_mom0, hi_proj = load_hi()
    co21, co21_proj = load_co21_IRAM()

    halpha_proj.quicklook()
    halpha_proj.FITSFigure.show_contour(hi_mom0, colors='w', alpha=0.8)
    halpha_proj.FITSFigure.show_contour(co21, cmap='autumn', alpha=0.8)

# Continuum, Halpha, CO21
elif whichplot == 5:

    continuum, continuum_proj = load_continuum_L()
    halp_kp, halp_kp_proj = load_halp_kp()
    co21, co21_proj = load_co21_IRAM()

    halp_kp_proj.quicklook()
    halp_kp_proj.FITSFigure.show_contour(continuum, colors='w', alpha=0.8)
    halp_kp_proj.FITSFigure.show_contour(co21, cmap='autumn', alpha=0.8)
