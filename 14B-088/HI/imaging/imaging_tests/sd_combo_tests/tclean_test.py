
'''
Given the same parameters, will tclean return the same result as clean?
'''

from tasks import tclean, clean, exportfits

# First make dirty images

tclean(vis='../14B-088_HI.ms.contsub',
       datacolumn='data',
       imagename='tclean_test_dirty_rightpointings_casa5.0',
       field='M33_2, M33_4, M33_1, M33_14',
       imsize=[2560, 2560],
       cell='3arcsec',
       specmode='cube',
       start=1001,
       width=1,
       nchan=1,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       restfreq='1420.40575177MHz',
       outframe='LSRK',
       pblimit=0.1,
       usemask='user',
       mask="M33_14B-088_HI_mask_channel_330.image",
       # pbmask=0.5,
       deconvolver='multiscale',
       scales=[]
      )


dirty_call_clean = clean(vis='../14B-088_HI.ms.contsub',
      imagename='clean_test_dirty_onlym33_14',
      # field='M33_1, M33_2, M33_4, M33_5, M33_6, M33_7_center, M33_8, M33_9, M33_10, M33_11, M33_12, M33_13, M33_14',
      field='M33_14'
      restfreq='1420.40575177MHz',
      mode='channel', width=1, nchan=1, start=1001,
      cell='3arcsec', multiscale=[],
      threshold='1.8mJy/beam', imagermode='mosaic', gain=0.1,
      imsize=[2560, 2560], weighting='natural', robust=0.0, niter=0,
      pbcor=False, minpb=0.1, interpolation='linear', usescratch=False,
      phasecenter='J2000 01h33m50.904 +30d39m35.79', veltype='radio',
      modelimage=None, mask="M33_14B-088_HI_mask_channel_330.image")

exportfits(imagename="tclean_test_dirty.image", fitsimage='tclean_test_dirty.image.fits',
           dropdeg=True)
exportfits(imagename="clean_test_dirty.image", fitsimage='clean_test_dirty.image.fits',
           dropdeg=True)

# Now clean them.

cleaned_call_tclean = tclean(vis='../14B-088_HI.ms.contsub',
       datacolumn='data',
       imagename='tclean_test_casa5.0',
       imsize=[2560, 2560],
       cell='3arcsec',
       specmode='cube',
       start=1001,
       width=1,
       nchan=1,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=100000,
       threshold='4mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       restfreq='1420.40575177MHz',
       outframe='LSRK',
       pblimit=0.1,
       usemask='user',
       mask="M33_14B-088_HI_mask_channel_330.image",
       # pbmask=0.5,
       deconvolver='multiscale',
       scales=[0, 6, 18]
      )


cleaned_call_clean = clean(vis='../14B-088_HI.ms.contsub',
      imagename='clean_test', field='M33*',
      restfreq='1420.40575177MHz',
      mode='channel', width=1, nchan=1, start=1001,
      cell='3arcsec', multiscale=[0, 6, 18],
      threshold='4mJy/beam', imagermode='mosaic', gain=0.1,
      imsize=[2560, 2560], weighting='natural', robust=0.0, niter=200000,
      pbcor=False, minpb=0.1, interpolation='linear', usescratch=False,
      phasecenter='J2000 01h33m50.904 +30d39m35.79', veltype='radio',
      modelimage=None, mask="M33_14B-088_HI_mask_channel_330.image")

exportfits(imagename="tclean_test.image", fitsimage='tclean_test.image.fits',
           dropdeg=True)
exportfits(imagename="clean_test.image", fitsimage='clean_test.image.fits',
           dropdeg=True)
exportfits(imagename="tclean_test.residual", fitsimage='tclean_test.residual.fits',
           dropdeg=True)
exportfits(imagename="clean_test.residual", fitsimage='clean_test.residual.fits',
           dropdeg=True)
exportfits(imagename="tclean_test.model", fitsimage='tclean_test.model.fits',
           dropdeg=True)
exportfits(imagename="clean_test.model", fitsimage='clean_test.model.fits',
           dropdeg=True)
