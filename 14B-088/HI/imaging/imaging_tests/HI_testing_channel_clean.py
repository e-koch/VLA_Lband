

import sys
import numpy as np

from tasks import clean

'''
Cleans an MS with a single channel given a mask and a model
'''

default("clean")

major, minor, revision = casadef.casa_version.split('.')
casa_version = 100 * int(major) + 10 * int(minor) + int(revision)

vis = sys.argv[-6]
model = sys.argv[-5]
mask = sys.argv[-4]
use_all_fields = True if sys.argv[-3] == "T" else False
use_multiscale = True if sys.argv[-2] == "T" else False
use_tclean = True if sys.argv[-1] == "T" else False

if model == "None":
    model = None
if mask == "None":
    mask = None

if use_tclean:
    if casa_version < 450:
        raise Warning("tclean only works for CASA versions above 4.5.")

if use_all_fields:
    field = 'M33*'
else:
    # Drop M33_3, the incorrect pointing.
    field = ",".join(["M33_{}".format(i) for i in range(1, 15)
                      if i not in [3, 7]]) + ", M33_7_center"

if use_multiscale:
    multiscale = [0, 4, 8, 20, 40, 80]
    # Different progression based on # pixels within a beam
    # multiscale = list(np.array([0, 1, 3, 9, 27, 81]) * 4)
    # multiscale = list(np.array([0, 2, 5]) * 4)
else:
    multiscale = []

out_root = "{0}.CASAVer_{1}.Model_{2}.Mask_{3}.AllFields_{4}.MScale_{5}" \
           ".Tclean_{6}".format(vis[:-3],
                                casa_version,
                                "T" if model is not None else "F",
                                "T" if mask is not None else "F",
                                "T" if use_all_fields else "F",
                                "T" if use_multiscale else "F",
                                "T" if use_tclean else "F")

if use_tclean:
    from tasks import tclean

    tclean(vis=vis, imagename=out_root + '.clean', field=field,
           restfreq='1420.40575177MHz', specmode='cube', nchan=1,
           start=1, width=1, cell='3arcsec', scales=multiscale,
           niter=200000, threshold="1.8mJy/bm", gain=0.1, imsize=[2560, 2560],
           gridder='mosaic', weighting='natural', veltype='radio', pblimit=0.2,
           interpolation='linear', startmodel=model, usemask='user', mask=mask,
           phasecenter='J2000 01h33m50.904 +30d39m35.79',
           )
else:
    clean(vis=vis, imagename=out_root + '.clean', field=field,
          restfreq='1420.40575177MHz',
          mode='channel', width=1, nchan=1, start=1,
          cell='3arcsec', multiscale=multiscale,
          threshold='1.8mJy/beam', imagermode='mosaic', gain=0.1,
          imsize=[2560, 2560], weighting='natural', robust=0.0, niter=200000,
          pbcor=True, minpb=0.2, interpolation='linear', usescratch=False,
          phasecenter='J2000 01h33m50.904 +30d39m35.79', veltype='radio',
          modelimage=model, mask=mask)
