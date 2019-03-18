
'''
Image per SPW the pipeline-calibrated data. Does statwt down-weight the
remaining outliers enough to avoid needing additional flagging?
'''

import sys
import os

from tasks import tclean

myvis = sys.argv[-1]

folder_name = 'mosaic_test_image'
if not os.path.exists(folder_name):
    os.mkdir(folder_name)

for spw_num in range(8):

    tclean(vis=myvis,
           datacolumn='corrected',
           imagename=os.path.join(folder_name,
                                  '{0}_spw_{1}.dirty'.format(myvis.rstrip(".ms"), spw_num)),
           field='M33*',
           spw="{}".format(spw_num),
           imsize=[3500, 3500],
           cell='2arcsec',
           specmode='mfs',
           startmodel=None,
           gridder='mosaic',
           weighting='natural',
           niter=0,
           threshold='0.1mJy/beam',
           phasecenter='J2000 01h33m50.904 +30d39m35.79',
           pblimit=0.2,
           usemask='pb',
           pbmask=0.2,
           deconvolver='hogbom',
           pbcor=False,
           interactive=False
           )

    tclean(vis=myvis,
           datacolumn='corrected',
           imagename=os.path.join(folder_name,
                                  '{0}_spw_{1}.clean_automask'.format(myvis.rstrip(".ms"), spw_num)),
           field='M33*',
           spw="{}".format(spw_num),
           imsize=[3500, 3500],
           cell='2arcsec',
           specmode='mfs',
           startmodel=None,
           gridder='mosaic',
           weighting='natural',
           niter=10000,
           threshold='0.1mJy/beam',
           phasecenter='J2000 01h33m50.904 +30d39m35.79',
           pblimit=0.2,
           usemask='auto-multithresh',
           pbmask=0.2,
           pbcor=False,
           interactive=False,
           deconvolver='multiscale',
           scales=[0, 6, 12],
           minpercentchange=2.,
           noisethreshold=4.,
           lownoisethreshold=2.,
           )
