
'''
Only image the region directly around the OH detection
'''

import os
import sys

from tasks import rmtables, impbcor, exportfits, tclean

# orig_dir = os.getcwd()

myspw = int(sys.argv[-4])

assert myspw in [3, 5, 6, 7]

use_contsub = True if sys.argv[-3] == "True" else False

myweighting = str(sys.argv[-2])
assert myweighting in ['natural', 'briggs', 'uniform']

# Run a single plane test to check spatial setup
single_chan_test = True if sys.argv[-1] == 'True' else False

# The VLA_Lband directory needs to be in the path
# from CASA_functions import (set_imagermode, set_imagesize, set_cellsize)

execfile(os.path.expanduser("~/ownCloud/code_development/VLA_Lband/17B-162/spw_setup.py"))


# With and without contsub
# for contsub in [True, False]:

spw = myspw

rest_freq = str(linespw_dict[spw][1])

line = linespw_dict[spw][0]

print("On line: {0} {1}".format(line, rest_freq))

# os.chdir(line)

if use_contsub:
    ms_active = "17B-162_{0}_spw_{1}_LSRK.ms.contsub".format(line, spw)
else:
    ms_active = "17B-162_{0}_spw_{1}_LSRK.ms".format(line, spw)

print("Running on {}".format(ms_active))

out_directory = "{}_imaging".format(line)

if not os.path.exists(out_directory):
    os.mkdir(out_directory)

default('tclean')

# Imaging parameters
# Common velocity range + width
imagermode = 'cube'
mode = 'velocity'
if single_chan_test:
    nchan = 1
else:
    # -80 to -300 km/s
    nchan = 167

width = "1.5km/s"
start = "-300.0km/s"
# threshold = "3mJy/beam"

# Fairly aggressive clean. Though this will be couple
# with automasking
nsigma = 2.0

# Imaging whole galaxy in B-pointings
field = "M33*"
# Location of known 1665 maser
# phasecenter = 'J2000 01h34m00.122 +30d40m47.19'
phasecenter = "J2000 01h33m50.904 +30d39m35.79"

spw = ""
imagermode = "mosaic"
outframe = "LSRK"
veltype = "radio"
pbcor = False
robust = 0.0
restfreq = rest_freq
interpolation = 'linear'
pblimit = 0.1

spw_num = 0
# max_size = 15000

# Determine imagermode, cell size, and image size
# imagermode = set_imagermode(ms_active, field.rstrip("*"))
# cellsize = set_cellsize(ms_active, 0, sample_factor=6.)
# imagesize = set_imagesize(ms_active, 0, field.rstrip("*"), sample_factor=6.,
#                           pblevel=pblimit, max_size=max_size)

if myweighting == 'natural':
    # Expect beams ~6"
    cellsize = "1.5arcsec"
    imagesize = 3200
# elif myweighting == 'briggs'
else:
    # Expect beam 3-4"
    cellsize = "0.75arcsec"
    imagesize = 6600

# cellsize = '0.65arcsec'
# imagesize = [7168, 6912]

print("Output image cell size: {}".format(cellsize))
print("Image pixel size: {}".format(imagesize))

if myweighting == 'briggs':
    weighting_label = 'robust_0'
else:
    weighting_label = myweighting

if single_chan_test:
    if use_contsub:
        myimgname = "{0}_17B-162_{1}_1_5kms_singlechantest".format(line, weighting_label)
    else:
        myimgname = "{0}_17B-162_wcont_{1}_1_5kms_singlechantest".format(line, weighting_label)
else:
    if use_contsub:
        myimgname = "{0}_17B-162_{1}_1_5kms".format(line, weighting_label)
    else:
        myimgname = "{0}_17B-162_wcont_{1}_1_5kms".format(line, weighting_label)

out_image = os.path.join(out_directory,
                         myimgname)

# Ensure there aren't any old versions
rmtables(out_image + ".*")

# Make a dirty map first
tclean(vis=ms_active,
       # datacolumn='data',
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
       weighting=myweighting,
       robust=robust,
       niter=0,
       # threshold=threshold,
       nsigma=nsigma,
       phasecenter=phasecenter,
       restfreq=rest_freq,
       outframe=outframe,
       pblimit=pblimit,
       usemask='auto-multithresh',
       sidelobethreshold=3.0,
       noisethreshold=5.0,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=True,  # All planes will be noise dominated.
       # mask='circle[[150pix,150pix],20pix]',
       deconvolver='hogbom',  # No extended structure expected.
       pbcor=False,
       chanchunks=-1,
       interactive=False,
       restoration=False,
       )

# Stop here for the single channel test
if single_chan_test:
    casalog.post("Only run single channel for single channel test.")
    sys.exit(0)

# Copy the dirty residual cube
orig_res_cube = "{}.residual".format(out_image)
backup_res_cube = "{}.residual".format(out_image)

os.system("cp -r {0} {1}".format(orig_res_cube, backup_res_cube))

tclean(vis=ms_active,
       # datacolumn='data',
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
       weighting=myweighting,
       robust=robust,
       niter=1000000,
       # threshold=threshold,
       nsigma=nsigma,
       phasecenter=phasecenter,
       restfreq=rest_freq,
       outframe=outframe,
       pblimit=pblimit,
       usemask='auto-multithresh',
       sidelobethreshold=3.0,
       noisethreshold=5.0,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=True,  # All planes will be noise dominated.
       # mask='circle[[150pix,150pix],20pix]',
       deconvolver='hogbom',  # No extended structure expected.
       pbcor=False,
       chanchunks=-1,
       interactive=False,
       restoration=True,
       calcres=False,
       calcpsf=False,
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
