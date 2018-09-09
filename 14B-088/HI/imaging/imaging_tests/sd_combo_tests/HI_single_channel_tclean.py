

import sys
import numpy as np
import os

from tasks import tclean, feather

'''
Cleans an MS with a single channel given a mask and a model
'''

default("tclean")

vis = sys.argv[-5]
model = sys.argv[-4]
use_model = True if sys.argv[-3] == "T" else False
mask = sys.argv[-2]
out_root = sys.argv[-1]

if model == "None":
    model = None
if mask == "None":
    mask = None

field = 'M33*'
multiscale = [0, 4, 8, 20, 40, 80]

if use_model:
    start_model = model

tclean(vis=vis, imagename=out_root + '.clean',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
       field=field, gridder='mosaic', weighting='natural',
       restfreq='1420.40575177MHz',
       specmode='cube', nchan=1, start=1, width=1,
       cell='3arcsec', imsize=[2560, 2560],
       deconvolver='multiscale', scales=multiscale,
       niter=200000, threshold="1.8mJy/bm", gain=0.1,
       veltype='radio', pblimit=0.1, pbcor=False,
       interpolation='linear', startmodel=start_model,
       usemask='user', mask=mask,
       )

# Run feathering with the model
if model is not None:
    if os.path.exists(out_root + ".clean.image"):
        feather(imagename=out_root + ".clean.image.feathered",
                highres=out_root + ".clean.image",
                lowres=model)
