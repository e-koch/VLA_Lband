
import sys

'''
Cleans an MS with a single channel given a mask and a model
'''

vis = sys.argv[4]
model = sys.argv[5]
mask = sys.argv[6]

out_root = vis[:-3]

clean(vis=vis, imagename=out_root+'.clean', field='M33*',
      restfreq='1420.40575177MHz',
      mode='channel', width=1, nchan=1, start=1,
      cell='1.5arcsec', multiscale=[0, 4, 8, 20, 40, 80],
      threshold='2.2mJy/beam', imagermode='mosaic', gain=0.2,
      imsize=[4096, 4096], weighting='natural', robust=0.0, niter=50000,
      pbcor=True, minpb=0.7, interpolation='linear', usescratch=True,
      cyclefactor=4,
      phasecenter='J2000 01h33m50.904 +30d39m35.79', veltype='radio',
      outframe='LSRK', modelimage=model, mask=mask)
