
'''
Parameters chosen based on the one OH(1665) detection.
'''

import os

from tasks import clean, rmtables, impbcor


# The VLA_Lband directory needs to be in the path
from CASA_functions import (set_imagermode, set_imagesize, set_cellsize)

default('clean')

ms_active = "OH1665_14B-088.ms.contsub"
out_directory = "imaging_1point5km_s"

# Imaging parameters
mode = 'velocity'
nchan = 15
width = "1.5km/s"
start = "-220.0km/s"
threshold = "1.0mJy"
field = "M33*"
phasecenter = 'J2000 01h33m50.904 +30d39m35.79'
spw = ""
imagermode = "mosaic"
multiscale = []
outframe = "LSRK"
veltype = "radio"
pbcor = False
robust = 0.0
restfreq = "1.6654018GHz"
interpolation = 'linear'
minpb = 0.1

spw_num = 0
max_size = 15000

# Determine imagermode, cell size, and image size
imagermode = set_imagermode(ms_active, field.rstrip("*"))
cellsize = set_cellsize(ms_active, 0, sample_factor=6.)
imagesize = set_imagesize(ms_active, 0, field.rstrip("*"), sample_factor=6.,
                          pblevel=minpb, max_size=max_size)

print("Output image cell size: {}".format(cellsize))
print("Image pixel size: {}".format(imagesize))

# Re-image with 3 different weighting choices
# robust has 0.0 set above
weightings = ["natural", "briggs", "uniform"]

for weighting in weightings:

    out_image = os.path.join(out_directory,
                             "OH1665_14B-088_{}".format(weighting))

    # Ensure there aren't any old versions
    rmtables(os.path.join(out_image, "*"))

    clean(vis=ms_active, imagename=out_image, niter=0, imsize=imagesize,
          cell=cellsize, mode=mode, nchan=nchan, width=width, start=start,
          threshold=threshold, field=field, phasecenter=phasecenter, spw=spw,
          imagermode=imagermode, multiscale=multiscale, outframe=outframe,
          veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
          restfreq=restfreq, usescratch=False,
          interpolation=interpolation, pbcor=True)

    # Now run the pbcor
    impbcor(imagename="{}.image".format(out_image),
            pbimage="{}.flux".format(out_image),
            outfile="{}.image.pbcor".format(out_image),
            cutoff=minpb, mode='divide')
