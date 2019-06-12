

import sys
import os
# from distutils.dir_util import mkpath
import re
import numpy as np
import time
import socket

from tasks import tclean, tget

'''
Cleans a single channel given the channel name
'''

# Load in the SPW dict in the repo on cedar
if socket.gethostname().lower() == 'segfault':
    execfile(os.path.expanduser("~/ownCloud/code_development/VLA_Lband/17B-162/spw_setup.py"))
else:
    execfile(os.path.expanduser("~/code/VLA_Lband/17B-162/spw_setup.py"))

chan_num = int(sys.argv[-3])

# Load in the imaging parameters from the given file name
parameter_file = sys.argv[-2]

# Assume file structure of channel_path/channel_${num}/
# Need to give the overall channel_path
channel_path = sys.argv[-1]

# Load parameters
tget(tclean, parameter_file)

# Append the full channel path to the vis's
# vis = [os.path.join(channel_path, "channel_{}".format(chan_num),
#                     "{0}_channel_{1}.ms".format(mss, chan_num))
#        for mss in vis]

vis = os.path.join(channel_path, "channel_{}".format(chan_num),
                   "14B_17B_channel_{}.ms".format(chan_num))

# Run initweights. Otherwise the C-config data is severely under-weighted
initweights(vis=vis, wtmode='nyq')

# Get the output path and create directories, if needed, based on
# the imagename
# if imagename.split("/") > 1:

#     output_path = "/".join(imagename.split("/")[:-1])

#     # Since the imaging is not run from the parent path for the channel output,
#     # use `mkpath`
#     if not os.path.exists(output_path):
#         mkpath(output_path)

# Now update the imagename with the channel number
imagename = os.path.join(channel_path, "channel_{}".format(chan_num),
                         "{0}_channel_{1}".format(imagename, chan_num))

# Based on the channel number, update the start velocity for this channel

# Get the value out from the string, removing the unit
try:
    split_start = filter(None, re.split(r'(\d+)', start))

    init_start = float("".join(split_start[:-1]))
    spec_unit = split_start[-1]

    chan_width = float("".join(filter(None, re.split(r'(\d+)', width))[:-1]))

    start_vel = init_start + chan_width * chan_num

    int_settings = False

except Exception:
    int_settings = True

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

    startmodel = os.path.join(channel_path,
                              "channel_".format(chan_num),
                              "{0}_channel_{1}.image"(startmodel.split(".image")[0],
                                                      chan_num))

    if not os.path.exists(startmodel):
        raise ValueError("Given startmodel does not exist")
else:
    startmodel = ""

# The naming scheme should split name.mask to name_channel_{}.mask
# The file MUST end in ".mask"
if mask is not None and len(mask) > 0 and usemask == "user":

    mask = os.path.join(channel_path,
                        "channel_{}".format(chan_num),
                        "{0}_channel_{1}".format(mask, chan_num))

    if not os.path.exists(mask):
        raise ValueError("Given mask name ({0}) does not exist".format(mask))

# Grab freq from the SPW dict
# spw_num = int(spw)
# We're only imaging HI here
spw_num = 0

# Only update a few parameters, as needed

if not int_settings:
    start = "{0}{1}".format(start_vel, spec_unit)
    width = "{0}{1}".format(chan_width, spec_unit)
    nchan = 1
else:
    start = 1
    width = 1
    nchan = 1

restfreq = linespw_dict[spw_num][1]
restart = True
calcres = do_calcres
calcpsf = do_calcpsf
interactive = 0  # Returns a summary dictionary

# out_dict = tclean()

## (1) Import the python application layer

# from imagerhelpers.imager_base import PySynthesisImager
# imagerInst = PySynthesisImager

from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
imagerInst = PyParallelContSynthesisImager

from imagerhelpers.input_parameters import ImagerParameters


# ## (2) Set up Input Parameters
# ## - List all parameters that you need here
# ## - Defaults will be assumed for unspecified parameters
# ## - Nearly all parameters are identical to that in the task. Please look at the
# ## list of parameters under __init__ using " help ImagerParameters " )


# Put all parameters into dictionaries and check them.
# Changes introduced in CASA 5.5
inpparams = locals().copy()
inpparams['msname'] = inpparams.pop('vis')
inpparams['timestr'] = inpparams.pop('timerange')
inpparams['uvdist'] = inpparams.pop('uvrange')
inpparams['obs'] = inpparams.pop('observation')
inpparams['state'] = inpparams.pop('intent')
inpparams['loopgain'] = inpparams.pop('gain')
inpparams['scalebias'] = inpparams.pop('smallscalebias')
defparm = dict(zip(ImagerParameters.__init__.__func__.__code__.co_varnames[1:],
                   ImagerParameters.__init__.func_defaults))
bparm = {k: inpparams[k] if k in inpparams else defparm[k] for k in defparm}
paramList = ImagerParameters(**bparm)

# (3) Construct the PySynthesisImager object, with all input parameters

# imager = PySynthesisImager(params=paramList)
imager = PyParallelContSynthesisImager(params=paramList)

# (4) Initialize various modules.
# - Pick only the modules you will need later on. For example, to only make
# the PSF, there is no need for the deconvolver or iteration control modules.

# Initialize modules major cycle modules
try:
    t0 = time.time()

    imager.initializeImagers()
    imager.initializeNormalizers()
    imager.setWeighting()

    t1 = time.time()

    casalog.post("Time for initializing imager and normalizers: " +
                 "%.2f" % (t1 - t0) + " sec")

    # Init minor cycle modules
    if niter > 0:
        t2 = time.time()
        imager.initializeDeconvolvers()
        t3 = time.time()
        casalog.post("Time for initializing deconvolver: " +
                     "%.2f" % (t3 - t2) + " sec")

    if niter > 0:
        t4 = time.time()
        imager.initializeIterationControl()
        t5 = time.time()
        casalog.post("Time for initializing iteration control: " +
                     "%.2f" % (t5 - t4) + " sec")

    # (5) Make the initial images

    if do_calcpsf:
        t6 = time.time()
        imager.makePSF()
        # From casa 5.5 tclean
        # if (psfphasecenter != '') and (gridder == 'mosaic'):
        #     # print "doing with different phasecenter psf"
        #     imager.unlockimages(0)
        #     psfParameters = paramList.getAllPars()
        #     psfParameters['phasecenter'] = psfphasecenter
        #     psfParamList = ImagerParameters(**psfParameters)
        #     psfimager = imagerInst(params=psfParamList)
        #     psfimager.initializeImagers()
        #     psfimager.setWeighting()
        #     psfimager.makeImage('psf', psfParameters['imagename']+'.psf')
        t7 = time.time()
        casalog.post("Time for creating PSF: " +
                     "%.2f" % (t7 - t6) + " sec")

        t8 = time.time()
        imager.makePB()
        t9 = time.time()
        casalog.post("Time for creating PB: " +
                     "%.2f" % (t9 - t8) + " sec")

    if do_calcres:
        casalog.post("Initial major cycle")

        t10 = time.time()
        imager.runMajorCycle()  # Make initial dirty / residual image
        t11 = time.time()
        casalog.post("Time for initial major cycle: " +
                     "%.2f" % (t11 - t10) + " sec")

        # Copy the initial residual map to a new name for post-imaging checks
        os.system("cp -r {0} {0}_init".format(imagename + ".residual"))

    if niter > 0:

        # (6) Make the initial clean mask
        imager.hasConverged()
        imager.updateMask()

        # (7) Run the iteration loops

        # Add an additional stopping criteria when the model flux between
        # major cycles changes by less than a set threshold.
        # Setting threshold to be 0.1%
        delta_model_flux_thresh = 1e-3

        model_flux_criterion = False

        mincyc_num = 0

        while not imager.hasConverged():
            # casalog.post("On minor cycle {}".format(mincyc_num))

            t0_l = time.time()
            imager.runMinorCycle()
            t1_l = time.time()
            casalog.post("Time for minor cycle: " +
                         "%.2f" % (t1_l - t0_l) + " sec")

            t2_l = time.time()
            imager.runMajorCycle()
            t3_l = time.time()
            casalog.post("Time for major cycle: " +
                         "%.2f" % (t3_l - t2_l) + " sec")

            summ = imager.IBtool.getiterationsummary()

            if mincyc_num == 0:
                model_flux_criterion = False
            else:
                model_flux_prev = summ['summaryminor'][2, :][-2]
                model_flux = summ['summaryminor'][2, :][-1]

                casalog.post("Previous model flux {0}. New model flux {1}".format(model_flux_prev, model_flux))

                model_flux_criterion = np.allclose(model_flux, model_flux_prev,
                                                   rtol=delta_model_flux_thresh)

            # Has the model converged?
            if model_flux_criterion:
                casalog.post("Model flux converged to within {}% between "
                             "major cycles.".format(delta_model_flux_thresh * 100))
                break
            else:
                time.sleep(10)

                imager.updateMask()

                time.sleep(10)

            mincyc_num += 1

    # (8) Finish up

    if niter > 0:
        out_dict = imager.IBtool.getiterationsummary()

        # Save the output dictionary. Numpy should be fine for this as the
        # individual channels will get concatenated together

        np.save(imagename + ".results_dict.npy", out_dict)

    if restoration:
        t12 = time.time()
        imager.restoreImages()
        t13 = time.time()
        casalog.post("Time for restoring images: " +
                     "%.2f" % (t13 - t12) + " sec")

        if pbcor:
            t14 = time.time()
            imager.pbcorImages()
            t15 = time.time()
            casalog.post("Time for pb-correcting images: " +
                         "%.2f" % (t15 - t14) + " sec")

    imager.deleteTools()

    t16 = time.time()
    casalog.post("Total Time: " +
                 "%.2f" % (t16 - t0) + " sec")


except Exception as e:
    casalog.post("Exception reported: {}".format(e), "SEVERE")
    casalog.post("Exception reported: {}".format(e.args), "SEVERE")

    try:
        imager.deleteTools()
    except Exception:
        pass

    raise e
