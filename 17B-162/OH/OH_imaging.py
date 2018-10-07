
'''
Only image the region directly around the OH detection
'''

import os

from tasks import rmtables, impbcor, exportfits, tclean

orig_dir = os.getcwd()


# The VLA_Lband directory needs to be in the path
# from CASA_functions import (set_imagermode, set_imagesize, set_cellsize)

# Loop through all 4 lines
spws = [3, 5, 6]#, 7]

execfile(os.path.expanduser("~/Dropbox/code_development/VLA_Lband/17B-162/spw_setup.py"))

for spw in spws:

    rest_freq = str(linespw_dict[spw][1])

    line = linespw_dict[spw][0]

    print("On line: {0} {1}".format(line, rest_freq))

    # os.chdir(line)

    ms_active = "17B-162_{0}_spw_{1}_LSRK.ms.contsub".format(line, spw)
    out_directory = "{}_imaging".format(line)

    if not os.path.exists(out_directory):
        os.mkdir(out_directory)

    # Imaging parameters
    mode = 'velocity'
    nchan = 45
    width = "1.5km/s"
    start = "-230.0km/s"
    # threshold = "3mJy/beam"
    nsigma = 2.0
    field = "M33*"
    phasecenter = 'J2000 01h34m00.122 +30d40m47.19'
    spw = ""
    imagermode = "mosaic"
    # multiscale = []
    outframe = "LSRK"
    veltype = "radio"
    pbcor = False
    robust = 0.0
    restfreq = rest_freq
    interpolation = 'linear'
    pblimit = 0.4

    spw_num = 0
    max_size = 15000

    # Determine imagermode, cell size, and image size
    # imagermode = set_imagermode(ms_active, field.rstrip("*"))
    # cellsize = set_cellsize(ms_active, 0, sample_factor=6.)
    # imagesize = set_imagesize(ms_active, 0, field.rstrip("*"), sample_factor=6.,
    #                           pblevel=pblimit, max_size=max_size)
    imagermode = 'cube'
    cellsize = '0.65arcsec'
    imagesize = [7168, 6912]

    print("Output image cell size: {}".format(cellsize))
    print("Image pixel size: {}".format(imagesize))

    # Re-image with 3 different weighting choices
    # robust has 0.0 set above
    # weightings = ["natural", "briggs", "uniform"]

    # weighting = "uniform"
    weighting = "natural"

    out_image = os.path.join(out_directory,
                             "{0}_17B-162_{1}".format(line, weighting))

    # Ensure there aren't any old versions
    rmtables(out_image + ".*")

    # This one is interactive so that a clean mask can be added.
    tclean(vis=ms_active,
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
           niter=10000,
           # threshold=threshold,
           nsigma=nsigma,
           phasecenter=phasecenter,
           restfreq=rest_freq,
           outframe=outframe,
           pblimit=pblimit,
           usemask='user',
           mask='circle[[150pix,150pix],20pix]',
           deconvolver='hogbom',
           pbcor=False,
           chanchunks=-1,
           interactive=False,
           )

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

    # os.chdir(orig_dir)
