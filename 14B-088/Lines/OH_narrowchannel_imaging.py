
'''
Parameters chosen based on the one OH(1665) detection.

To be run in mpicasa
'''

import os
import json
import time

from tasks import rmtables, impbcor, exportfits
from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
from imagerhelpers.input_parameters import ImagerParameters


orig_dir = os.getcwd()


# The VLA_Lband directory needs to be in the path
# from CASA_functions import (set_imagermode, set_imagesize, set_cellsize)

# Loop through all 4 lines

oh_lines = ['OH1612', 'OH1665', 'OH1667', 'OH1720']

line_file = \
    os.path.expanduser(
        "~/Dropbox/code_development/VLA_Lband/14B-088/Lines/spw_dict.txt")
with open(line_file, 'r') as f:
    line_dict = json.load(f)

for line in oh_lines:
    rest_freq = str(line_dict[line])

    print("On line: {0} {1}".format(line, rest_freq))

    os.chdir(line)

    ms_active = "{0}_14B-088.ms.contsub".format(line)
    out_directory = "imaging_1point5km_s".format(line)

    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    # Imaging parameters
    mode = 'velocity'
    nchan = 135
    width = "1.5km/s"
    start = "-280.0km/s"
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
    restfreq = rest_freq
    interpolation = 'linear'
    pblimit = 0.1

    spw_num = 0
    max_size = 15000

    # Determine imagermode, cell size, and image size
    # imagermode = set_imagermode(ms_active, field.rstrip("*"))
    # cellsize = set_cellsize(ms_active, 0, sample_factor=6.)
    # imagesize = set_imagesize(ms_active, 0, field.rstrip("*"), sample_factor=6.,
    #                           pblevel=pblimit, max_size=max_size)
    imagermode = 'cube'
    cellsize = '2.3arcsec'
    imagesize = [2560, 3000]

    print("Output image cell size: {}".format(cellsize))
    print("Image pixel size: {}".format(imagesize))

    # Re-image with 3 different weighting choices
    # robust has 0.0 set above
    # weightings = ["natural", "briggs", "uniform"]

    weightings = ["natural", "uniform"]

    for weighting in weightings:

        out_image = os.path.join(out_directory,
                                 "{0}_14B-088_{1}".format(line, weighting))

        # Ensure there aren't any old versions
        rmtables(out_image + ".*")

        paramList = \
            ImagerParameters(msname=ms_active,
                             datacolumn='data',
                             field=field,
                             imagename=out_image,
                             imsize=imagesize,
                             cell=cellsize,
                             specmode='cube',
                             start=start,
                             width=width,
                             nchan=nchan,
                             startmodel=None,
                             gridder='mosaic',
                             weighting=weighting,
                             niter=0,
                             threshold=threshold,
                             phasecenter=phasecenter,
                             restfreq=rest_freq,
                             outframe=outframe,
                             pblimit=pblimit,
                             usemask='pb',
                             mask=None,
                             deconvolver='hogbom',
                             dopbcorr=False,
                             chanchunks=-1
                             )

        imager = PyParallelCubeSynthesisImager(params=paramList)

        # init major cycle elements
        imager.initializeImagers()
        imager.initializeNormalizers()
        imager.setWeighting()

        # Init minor cycle elements
        imager.initializeDeconvolvers()
        imager.initializeIterationControl()
        imager.makePSF()
        imager.makePB()

        # Make dirty image
        t0 = time.time()
        imager.runMajorCycle()
        t1 = time.time()
        casalog.post("Time for major cycle: {}".format(t1 - t0))

        # retrec = imager.getSummary()

        imager.restoreImages()
        # imager.pbcorImages()
        concattype = 'virtualcopy'
        imager.concatImages(type=concattype)
        imager.deleteTools()

        # Now run the pbcor
        impbcor(imagename="{}.image".format(out_image),
                pbimage="{}.pb".format(out_image),
                outfile="{}.image.pbcor".format(out_image),
                cutoff=pblimit, mode='divide')

        exportfits(imagename="{}.image.pbcor".format(out_image),
                   fitsimage="{}.image.pbcor.fits".format(out_image),
                   velocity=True, overwrite=True,
                   dropdeg=True, history=False)

        exportfits(imagename="{}.pb".format(out_image),
                   fitsimage="{}.pb.fits".format(out_image),
                   velocity=True, overwrite=True,
                   dropdeg=True, history=False)

    os.chdir(orig_dir)
