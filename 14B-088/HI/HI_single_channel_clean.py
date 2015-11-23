

import sys

from casa_tools import myclean

'''
Cleans an MS with a single channel given a mask and a model
'''

vis = sys.argv[-3]
model = sys.argv[-2]
mask = sys.argv[-1]

if model == "None":
    model = None
if mask == "None":
    mask = None

out_root = vis[:-3]

myclean(vis=vis, imagename=out_root+'.clean', field='M33*',
        restfreq='1420.40575177MHz',
        mode='channel', width=1, nchan=1, start=1,
        cell='3arcsec', multiscale=[0, 4, 8, 20, 40, 80],
        threshold='1.5mJy/beam', imagermode='mosaic', gain=0.1,
        imsize=[2560, 2560], weighting='natural', robust=0.0, niter=200000,
        pbcor=True, minpb=0.2, interpolation='linear', usescratch=False,
        phasecenter='J2000 01h33m50.904 +30d39m35.79', veltype='radio',
        modelimage=model, mask=mask)
