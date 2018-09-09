
'''
Imaging tests for the 14B-088 continuum (I) data.
'''

import os

from tasks import tclean

vis = "14B-088_continuum_I.ms"
output_path = "imaging_nosub"

if not os.path.exists(output_path):
    os.mkdir(output_path)


tclean(vis=vis,
       datacolumn='data',
       imagename=os.path.join(output_path, 'M33_14B-088_continuum.dirty'),
       field='M33*',
       spw="1",
       imsize=[2560, 2560],
       cell='3arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       pblimit=-1,
       usemask='pb',
       pbmask=0.2,
       deconvolver='hogbom',
       pbcor=False,
       interactive=True
       )
