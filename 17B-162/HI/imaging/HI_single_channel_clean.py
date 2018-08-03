

import sys
import os
from distutils.dir_util import mkpath
import re
import numpy as np

from tasks import tclean, tget

'''
Cleans a single channel given the channel name
'''

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/17B-162/spw_setup.py"))
# execfile(os.path.expanduser("~/Dropbox/code_development/VLA_Lband/17B-162/spw_setup.py"))

chan_num = int(sys.argv[-2])

# Load in the imaging parameters from the given file name
parameter_file = sys.argv[-1]

# Load parameters
tget(tclean, parameter_file)

# Get the output path and create directories, if needed, based on
# the imagename
if imagename.split("/") > 1:

    output_path = "/".join(imagename.split("/")[:-1])

    # Since the imaging is not run from the parent path for the channel output,
    # use `mkpath`
    if not os.path.exists(output_path):
        mkpath(output_path)

# Now update the imagename with the channel number
imagename = "{0}_channel_{1}".format(imagename, chan_num)

# Based on the channel number, update the start velocity for this channel

# Get the value out from the string, removing the unit
split_start = filter(None, re.split(r'(\d+)', start))

init_start = float("".join(split_start[:-1]))
spec_unit = split_start[-1]

chan_width = float("".join(filter(None, re.split(r'(\d+)', width))[:-1]))

start_vel = init_start + chan_width * chan_num

# Check if the products already exist and we should avoid recomputing the PSF
# and residuals

if os.path.exists("{}.image".format(imagename)):
    casalog.post("Image already exists! Assuming this channel is completed.")
    import sys
    sys.exit(0)

# Check for the PSF and assume the rest of the products are there too.
if os.path.exists("{}.psf".format(imagename)):
    do_calcres = False
    do_calcpsf = False

    # Force startmodel to use the model on disk
    startmodel = None

else:
    do_calcres = True
    do_calcpsf = True

# If model or mask names are given, ensure they exist.
# These should already be split into individual channels for use here
# The naming scheme should split imagename.image to imagename_channel_{}.image
# The file MUST end in ".image"
if startmodel is not None and len(startmodel) > 0:

    startmodel = "{0}_channel_{1}.image"(startmodel.split(".image")[0],
                                         chan_num)

    if not os.path.exists(startmodel):
        raise ValueError("Given startmodel does not exist")

# The naming scheme should split name.mask to name_channel_{}.mask
# The file MUST end in ".mask"
if mask is not None and len(mask) > 0 and usemask == "user":

    mask = "{0}_channel_{1}.mask"(mask.split(".image")[0], chan_num)

    if not os.path.exists(mask):
        raise ValueError("Given mask name ({0}) does not exist".format(mask))

# Grab freq from the SPW dict
spw_num = int(spw)

# Only update a few parameters, as needed

start = "{0}{1}".format(start_vel, spec_unit)
width = "{0}{1}".format(chan_width, spec_unit)
nchan = 1
restfreq = linespw_dict[spw_num][1]
restart = True
calcres = do_calcres
calcpsf = do_calcpsf
interactive = 0  # Returns a summary dictionary

# out_dict = tclean()

# Alternatively, skip the mandatory plotting stage from the summary call

## (1) Import the python application layer

from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

## (2) Set up Input Parameters
## - List all parameters that you need here
## - Defaults will be assumed for unspecified parameters
## - Nearly all parameters are identical to that in the task. Please look at the
## list of parameters under __init__ using " help ImagerParameters " )


    # Put all parameters into dictionaries and check them.
paramList = ImagerParameters(msname=vis,
                             field=field,
                             spw=spw,
                             timestr=timerange,
                             uvdist=uvrange,
                             antenna=antenna,
                             scan=scan,
                             obs=observation,
                             state=intent,
                             datacolumn=datacolumn,
                             imagename=imagename,
                             imsize=imsize,
                             cell=cell,
                             phasecenter=phasecenter,
                             stokes=stokes,
                             projection=projection,
                             startmodel=startmodel,
                             specmode=specmode,
                             reffreq=reffreq,
                             nchan=1,
                             start="{0}{1}".format(start_vel, spec_unit),
                             width="{0}{1}".format(chan_width, spec_unit),
                             outframe=outframe,
                             veltype=veltype,
                             restfreq=linespw_dict[spw_num][1],
                             sysvel='',
                             sysvelframe='',
                             interpolation=interpolation,
                             gridder=gridder,
                             facets=facets,
                             chanchunks=chanchunks,
                             wprojplanes=wprojplanes,
                             vptable=vptable,
                             aterm=aterm,
                             psterm=psterm,
                             wbawp = wbawp,
                             cfcache = cfcache,
                             conjbeams = conjbeams,
                             computepastep =computepastep,
                             rotatepastep = rotatepastep,
                             pblimit=pblimit,
                             normtype=normtype,
                             outlierfile=outlierfile,
                             restart=restart,
                             weighting=weighting,
                             robust=robust,
                             npixels=npixels,
                             uvtaper=uvtaper,
                             niter=niter,
                             cycleniter=cycleniter,
                             loopgain=gain,
                             threshold=threshold,
                             nsigma=0.0,
                             cyclefactor=cyclefactor,
                             minpsffraction=minpsffraction,
                             maxpsffraction=maxpsffraction,
                             interactive=0,
                             deconvolver=deconvolver,
                             scales=scales,
                             nterms=nterms,
                             scalebias=smallscalebias,
                             restoringbeam=restoringbeam,
                             usemask=usemask,
                             mask=mask,
                             pbmask=pbmask,
                             sidelobethreshold=sidelobethreshold,
                             noisethreshold=noisethreshold,
                             lownoisethreshold=lownoisethreshold,
                             negativethreshold=negativethreshold,
                             smoothfactor=smoothfactor,
                             minbeamfrac=minbeamfrac,
                             cutthreshold=cutthreshold,
                             growiterations=growiterations,
                             dogrowprune=True,
                             minpercentchange=0.0,
                             verbose=True,
                             savemodel=savemodel,
                             )

# (3) Construct the PySynthesisImager object, with all input parameters

imager = PySynthesisImager(params=paramList)

# (4) Initialize various modules.
# - Pick only the modules you will need later on. For example, to only make
# the PSF, there is no need for the deconvolver or iteration control modules.

# Initialize modules major cycle modules

imager.initializeImagers()
imager.initializeNormalizers()
imager.setWeighting()

# Init minor cycle modules

if restoration and niter > 0:
    imager.initializeDeconvolvers()

if niter > 0:
    imager.initializeIterationControl()

# (5) Make the initial images

if do_calcpsf:
    imager.makePSF()

    imager.makePB()

if do_calcres:
    imager.runMajorCycle()  # Make initial dirty / residual image

if niter > 0:
    # (6) Make the initial clean mask
    imager.hasConverged()
    imager.updateMask()

    # (7) Run the iteration loops
    mincyc_num = 1
    while (not imager.hasConverged()):
        casalog.post("On minor cycle {}".format(mincyc_num))
        imager.runMinorCycle()
        imager.runMajorCycle()
        imager.updateMask()
        mincyc_num += 1

# (8) Finish up

if niter > 0:
    out_dict = imager.IBtool.getiterationsummary()

    # Save the output dictionary. Numpy should be fine for this as the
    # individual channels will get concatenated together

    np.save(imagename + ".results_dict.npy", out_dict)

if restoration:
    imager.restoreImages()

    if pbcor:
        imager.pbcorImages()

imager.deleteTools()
