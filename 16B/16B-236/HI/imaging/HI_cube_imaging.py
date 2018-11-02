
'''
Image only the A-array HI data with moderate cleaning (to start)
'''

import os

from tasks import tclean

execfile(os.path.expanduser("~/code/VLA_Lband/16B/spw_setup.py"))

myvis = '16B-236_HI_spw_0_LSRK.ms.contsub'

spw_num = 0

imagename = 'M33_16B-236_{0}_spw_{1}.clean_5sig'\
    .format(linespw_dict[spw_num][0], spw_num)

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
       spw=str(spw_num),
       field='M33*',
       imsize=[8000, 8000],
       cell="0.3arcsec",
       specmode='cube',
       start="-270km/s",
       width="5km/s",
       nchan=40,
       startmodel=None,
       gridder='standard',
       weighting='natural',
       niter=1000000,
       nsigma=5,
       phasecenter='J2000 01h33m33.191 +30d32m06.72',
       restfreq=linespw_dict[spw_num][1],
       outframe='LSRK',
       pblimit=0.5,
       usemask='pb',
       mask=None,
       deconvolver='multiscale',
       scales=[0, 6, 18],
       pbcor=False,
       veltype='radio',
       chanchunks=-1,
       restoration=True,
       parallel=True,
       calcpsf=calcpsf,
       calcres=calcres,
       )
