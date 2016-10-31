

import sys
import os

# from casa_tools import myclean
from tasks import clean, feather, rmtables

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

# Define a valid mask. This stays true so long as a given mask contains a
# usable region.
valid_mask = True
if mask is not None:
    # Check if there's anything in the mask before cleaning
    ia.open(mask)
    stats_dict = ia.statistics()
    ia.close()
    # If there's nothing there, max == min
    max_val = stats_dict["max"]
    if max_val == 0:
        casalog.post("No valid region in the given mask. Skipping clean.")
        valid_mask = False
    else:
        casalog.post("Mask contains valid region. Proceeding with clean.")

if valid_mask:
    clean(vis=vis, imagename=out_root + '.clean', field='M33*',
          restfreq='1420.40575177MHz',
          mode='channel', width=1, nchan=1, start=1,
          cell='3arcsec', multiscale=[0, 4, 8, 20, 40, 80],
          threshold='4.5mJy/beam', imagermode='mosaic', gain=0.1,
          imsize=[2560, 2560], weighting='briggs', robust=0.0, niter=200000,
          pbcor=True, minpb=0.2, interpolation='linear', usescratch=False,
          phasecenter='J2000 01h33m50.904 +30d39m35.79', veltype='radio',
          modelimage=model, mask=mask)

    # Only run if the image was actually produces
    if os.path.exists(out_root + ".clean.image"):
        feather(imagename=out_root + ".clean.image.feathered",
                highres=out_root + ".clean.image",
                lowres=model)

    # If something went awry, and the image wasn't produced, remove the
    # remnants.
    if not os.path.exists(out_root + ".clean.image"):
        rmtables(out_root + ".clean.*")
