

import sys

from casa_tools import myclean, myconcat

'''
Cleans an MS with a single channel given a mask and a model
'''

vis_1 = sys.argv[-5]
vis_2 = sys.argv[-4]
model = sys.argv[-3]
mask = sys.argv[-2]
out_root = sys.argv[-1]

if model == "None":
    model = None
if mask == "None":
    mask = None

myconcat(vis=[vis_1, vis_2], output_vis=out_root+".ms")

myclean(vis=out_root+".ms", imagename=out_root+'.clean', field='M33*',
        restfreq='1420.40575177MHz',
        mode='channel', width=1, nchan=1, start=1,
        cell='1.5arcsec', multiscale=[0, 4, 8, 20, 40, 80, 160],
        threshold='3.5mJy/beam', imagermode='mosaic', gain=0.1,
        imsize=[4096, 4096], weighting='briggs', robust=0.0, niter=200000,
        pbcor=True, minpb=0.2, interpolation='linear', usescratch=False,
        phasecenter='J2000 01h33m50.904 +30d39m35.79', veltype='radio',
        modelimage=model, mask=mask)
