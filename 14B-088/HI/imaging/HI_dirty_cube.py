
'''
Create a dirty HI cube for comparison and use in feathering.

*Note*: Ran with casa-prerelease-5.0.0-187.el7 to take advantage of tclean's
read-only mode, which speeds things up considerably.
'''

import os

from tasks import tclean, impbcor

# CASA init should have the VLA_Lband repo appended to the path
from paths import data_path

# The full MS isn't in the same path as the usual data products.
# full_path = os.path.join(data_path, "../VLA/14B-088/Lines/HI/")

# output_path = os.path.join(full_path, "dirty_cube")
output_path = "dirty_cube"


if not os.path.exists(output_path):
    os.mkdir(output_path)


# Image ALL channels in the continuum subtracted MS (~2000).
# Keep the same spatial settings as is used for the cleaned cubes.

# tclean(vis=os.path.join(full_path, '14B-088_HI.ms.contsub'),
tclean(vis='14B-088_HI_LSRK.ms.contsub',
       datacolumn='data',
       imagename=os.path.join(output_path, 'M33_14B-088_HI.contsub.dirty.600_1000'),
       field='M33*',
       imsize=[2560, 2560],
       cell='3arcsec',
       specmode='cube',
       start=600,
       width=1,
       nchan=400,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       restfreq='1420.40575177MHz',
       outframe='LSRK',
       pblimit=0.05,
       usemask='pb',
       pbmask=0.1,
       deconvolver='hogbom',
       pbcor=False,
       chanchunks=-1
       )

tclean(vis='14B-088_HI_LSRK.ms.contsub',
       datacolumn='data',
       imagename=os.path.join(output_path, 'M33_14B-088_HI.contsub.dirty.1001_1400'),
       field='M33*',
       imsize=[2560, 2560],
       cell='3arcsec',
       specmode='cube',
       start=1001,
       width=1,
       nchan=400,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       restfreq='1420.40575177MHz',
       outframe='LSRK',
       pblimit=0.05,
       usemask='pb',
       pbmask=0.1,
       deconvolver='hogbom',
       pbcor=False,
       chanchunks=-1
       )

tclean(vis='14B-088_HI_LSRK.ms.contsub',
       datacolumn='data',
       imagename=os.path.join(output_path, 'M33_14B-088_HI.contsub.dirty.1401_1800'),
       field='M33*',
       imsize=[2560, 2560],
       cell='3arcsec',
       specmode='cube',
       start=1401,
       width=1,
       nchan=400,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       restfreq='1420.40575177MHz',
       outframe='LSRK',
       pblimit=0.05,
       usemask='pb',
       pbmask=0.1,
       deconvolver='hogbom',
       pbcor=False,
       chanchunks=-1
       )

tclean(vis='14B-088_HI_LSRK.ms.contsub',
       datacolumn='data',
       imagename=os.path.join(output_path, 'M33_14B-088_HI.contsub.dirty.1801_2001'),
       field='M33*',
       imsize=[2560, 2560],
       cell='3arcsec',
       specmode='cube',
       start=1801,
       width=1,
       nchan=-1,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       restfreq='1420.40575177MHz',
       outframe='LSRK',
       pblimit=0.05,
       usemask='pb',
       pbmask=0.1,
       deconvolver='hogbom',
       pbcor=False,
       chanchunks=-1
       )


# Apply pb correction
# impbcor(imagename=os.path.join(output_path, 'M33_14B-088_HI.dirty.image'),
        # pbimage=os.path.join(output_path, 'M33_14B-088_HI.dirty.pb'),
        # outfile=os.path.join(output_path, 'M33_14B-088_HI.dirty.image.pbcor'))
