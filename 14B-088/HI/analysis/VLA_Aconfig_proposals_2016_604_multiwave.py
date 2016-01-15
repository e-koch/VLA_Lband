
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

hi_mom0 = \
    fits.open("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"
              "M33_14B-088_HI.clean.image.mom0_3sigmask.fits")[0]

hi_proj = Projection(hi_mom0.data, unit=u.Jy, wcs=WCS(hi_mom0.header))

halpha = fits.open("/media/eric/Data_3/M33/Hubble/hst_05773_03_wfpc2_f656n_wf"
                   "/hst_05773_03_wfpc2_f656n_wf_sci.fits")[1]

halpha_proj = Projection(halpha.data, unit=u.count, wcs=WCS(halpha.header))


sii = fits.open("/media/eric/Data_3/M33/Hubble/hst_05773_03_wfpc2_f673n_wf/"
                "hst_05773_03_wfpc2_f673n_wf_sci.fits")[1]

sii_proj = Projection(sii.data, unit=u.count, wcs=WCS(sii.header))

chandra = fits.open("/media/eric/Data_3/M33/Chandra/images/"
                    "broad_thresh_exp_sm.img")[0]

chandra_proj = Projection(chandra.data, unit=u.count, wcs=WCS(chandra.header))

co21 = fits.open("/media/eric/Data_3/M33/IRAM/CO_12sec_mom0.fits")[0]

new_wcs = drop_axis(drop_axis(WCS(co21.header), 3), 2)
co21_proj = Projection(co21.data.squeeze(), unit=u.K, wcs=new_wcs)


halp_sd = fits.open("/media/eric/Data_3/M33/Halpha/ha.fits")[0]

halp_sd_proj = Projection(halp_sd.data, unit=u.count, wcs=WCS(halp_sd.header))

continuum = fits.open("/media/eric/Data_2/14B-088_continuum/test_imaging/"
                      "NGC_604_M33_3_briggs_0.5.image.tt0.fits")[0]

new_wcs = drop_axis(drop_axis(WCS(continuum.header), 3), 2)

continuum_proj = Projection(continuum.data.squeeze(), unit=u.Jy,
                            wcs=new_wcs)


# HI, Halpha, 2-10 keV
if whichplot == 0:

    halpha_proj.quicklook()
    halpha_proj.FITSFigure.show_contour(hi_mom0, colors='w', alpha=0.8)
    halpha_proj.FITSFigure.show_contour(chandra, colors='r', alpha=0.8)

# HI, [SII], 2-10 keV
elif whichplot == 1:

    sii_proj.quicklook()
    sii_proj.FITSFigure.show_contour(hi_mom0, colors='w', alpha=0.8)
    sii_proj.FITSFigure.show_contour(chandra, colors='r', alpha=0.8, )

# HI, Halpha SD, CO21
elif whichplot == 2:

    halp_sd_proj.quicklook()
    halp_sd_proj.FITSFigure.show_contour(hi_mom0, colors='g', alpha=0.8)
    halp_sd_proj.FITSFigure.show_contour(co21, colors='r', alpha=0.8)

# HI, Halpha, CO21
elif whichplot == 3:

    halpha_proj.quicklook()
    halpha_proj.FITSFigure.show_contour(hi_mom0, colors='w', alpha=0.8)
    halpha_proj.FITSFigure.show_contour(co21, cmap='autumn', alpha=0.8)

# Continuum, Halpha, CO21
elif whichplot == 4:

    halp_sd_proj.quicklook()
    halp_sd_proj.FITSFigure.show_contour(continuum, colors='w', alpha=0.8)
    halp_sd_proj.FITSFigure.show_contour(co21, cmap='autumn', alpha=0.8)
