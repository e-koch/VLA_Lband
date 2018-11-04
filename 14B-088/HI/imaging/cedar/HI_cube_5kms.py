
'''
Make a wide channel 14B cube.
'''

import os

from tasks import tclean


myvis = '14B-088_HI_LSRK.ms.contsub'

imagename = 'M33_14B-088_HI_5kms.clean'

if os.path.exists("{}.psf".format(imagename)):
    calcres = False
    calcpsf = False

    # Force startmodel to use the model on disk
    startmodel = None

else:
    calcres = True
    calcpsf = True

tclean(vis=myvis,
       datacolumn='corrected',
       imagename=imagename,
       spw='0',
       field='M33*',
       imsize=[2560, 2560],
       cell="3arcsec",
       specmode='cube',
       start="-330km/s",
       width="5km/s",
       nchan=58,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=1000000,
       nsigma=2.,
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       restfreq='1420.40575177MHz',
       outframe='LSRK',
       pblimit=0.3,
       usemask='pb',
       mask=None,
       deconvolver='multiscale',
       scales=[0, 6, 12, 30, 60, 120],
       pbcor=False,
       veltype='radio',
       chanchunks=-1,
       restoration=True,
       parallel=True,
       calcpsf=calcpsf,
       calcres=calcres,
       )
