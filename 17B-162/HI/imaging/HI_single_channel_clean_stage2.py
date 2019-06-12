

import sys
import os
import re
import numpy as np
import time
import scipy.ndimage as nd
from glob import glob
import socket

from tasks import tclean, tget

'''
Stage 2 of cleaning single channels.

Expects that the data have been cleaned to 5-sigma with
HI_single_channel_clean.py

Stage 2 creates a clean mask based on the 5-sigma-cleaned image by identifying
the emission footprint. The cleaning is continued on the masked region to
~2-sigma.

It's likely best to limit the largest scales used for multi-scale clean here,
unlike the first stage.
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
#                     "{0}_channel_{1}".format(mss, chan_num))
#        for mss in vis]

vis = os.path.join(channel_path, "channel_{}".format(chan_num),
                   "14B_17B_channel_{}.ms".format(chan_num))

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

# Check if the summary dictionary has already been saved for this channel
# If so, there is no need to clean any further

summ_files = glob("{}*.npy".format(imagename))

if len(summ_files) == 2:
    casalog.post("Summary file already exists! No need to clean further!")
    import sys
    sys.exit(0)

do_calcres = False
do_calcpsf = False

# Force startmodel to use the model on disk
startmodel = None


# Open the image, threshold at ~2-sigma, reject small regions,
# then expand the remaining regions to ensure all emission is captured in the
# masked region.

# myimage = imagename + ".image"

# new_maskname = "{0}.mask_stage2".format(imagename)

# if not os.path.exists(new_maskname):

#     # Also set the name of the mask
#     maskname = "{0}.mask".format(imagename)

#     ia.open(maskname)

#     orig_mask_data = ia.getchunk().squeeze()

#     ia.done()
#     ia.close()

#     # arcsec
#     major = imhead(imagename=myimage, mode='get', hdkey='beammajor')['value']
#     minor = imhead(imagename=myimage, mode='get', hdkey='beamminor')['value']
#     unitCDELT = imhead(imagename=myimage, mode='get', hdkey='cdelt2')['unit']
#     cdelt2 = imhead(imagename=myimage, mode='get', hdkey='cdelt2')['value']
#     assert unitCDELT == 'rad'
#     pixelsize = cdelt2 * (180. / np.pi) * 3600
#     beamarea = (major * minor * np.pi) / (pixelsize**2)

#     # Create new thresholded mask

#     # Threshold set by comparing to the cycle threshold saved in the previous
#     # clean call
#     results_dict = np.load(imagename + ".results_dict.npy").item()

#     # Set mask threshold at 2 * sigma
#     # Grab the residual peak after the last minor cycle. This
#     # is assumed to be set to 5-sigma from the first cycle!
#     sigma = (results_dict['summaryminor'][1][-1] / 5.)
#     thresh = sigma

#     immath(imagename=myimage, outfile=new_maskname,
#            expr='iif(IM0 > ' + str(thresh) + ', 1.0, 0.0)')

#     casalog.post("Creating clean mask")
#     casalog.post("Beam area is {0:.2f} pix^2".format(beamarea))

#     # Make a second mask, that will not be saved, at a threshold
#     # twice the original. Good regions must contain one point above
#     # the higher threshold
#     ia.open(myimage)

#     high_thresh_mask = ia.getchunk().squeeze() > thresh * 4.

#     ia.done()
#     ia.close()

#     # Remove small regions
#     ia.open(new_maskname)
#     mask_data = ia.getchunk()

#     orig_shape = mask_data.shape

#     mask_data = mask_data.squeeze()

#     labels_mask, num = nd.label(mask_data, np.ones((3, 3)))

#     pixel_sizes = nd.sum(mask_data, labels_mask, range(1, num + 1))

#     casalog.post("Removing small regions from clean mask")
#     casalog.post("Removing {0} of {1} regions from mask"
#                  .format(sum(pixel_sizes < beamarea * 0.5), num + 1))

#     new_mask_data = np.zeros_like(mask_data, dtype=bool)

#     good_labels = np.where(pixel_sizes > beamarea * 0.5)[0] + 1

#     for label in good_labels:

#         posns = np.where(labels_mask == label)

#         # Check if any points in the region reach the higher
#         # threshold. If not, don't include in the mask
#         # Set the minimum at a tenth of the beam. Avoid single
#         # high-noise spikes
#         if sum(high_thresh_mask[posns]) < 0.3 * beamarea:
#             continue

#         new_mask_data[posns] = True

#     mask_data = new_mask_data

#     # Check whether anything is left. If there isn't, there's no need
#     # to clean anymore!
#     if not mask_data.any():
#         casalog.post("No regions contained in the mask! No further cleaning"
#                      " is required.")


#         ia.putchunk(mask_data.reshape(orig_shape).astype(int))
#         ia.done()
#         ia.close()

#         import sys
#         sys.exit(0)

#     casalog.post("Smoothing clean mask by six times larger than largest scale.")

#     # Grab the largest multiscale scale
#     dilate_size = scales[-1] * 6
#     radius = dilate_size / 2
#     L = np.arange(-radius, radius + 1)
#     X, Y = np.meshgrid(L, L)
#     dilate_element = np.array((X ** 2 + Y ** 2) <= radius ** 2)

#     mask_data = nd.binary_dilation(mask_data, structure=dilate_element,
#                                    mask=orig_mask_data.astype(bool),
#                                    iterations=5)

#     # We should largely expect a single massive blob in the mask since the
#     # shape of the emission is dominated by the circular rotation
#     # So it is fairly safe to fill in holes in the mask
#     mask_data = nd.binary_fill_holes(mask_data)

#     ia.putchunk(mask_data.reshape(orig_shape).astype(int))
#     ia.done()
#     ia.close()

#     casalog.post("Finished creating clean mask")

# else:
#     casalog.post("Found a clean mask. Skipping making a new one.")

# mask = new_maskname
# usemask = "user"

# Grab freq from the SPW dict
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

# print(argh)

# Copy the first stage files into another directory, in case the deep
# clean goes badly
keep_old_path = os.path.join(channel_path, "channel_{}".format(chan_num),
                             "stage_1")
if os.path.exists(keep_old_path):
    os.system("rm -r {}".format(keep_old_path))

os.mkdir(keep_old_path)

stage1_files = glob(imagename + "*")

for fi in stage1_files:
    os.system("cp -r {0} {1}".format(fi, keep_old_path))


# Don't do this if we're using a signal mask from the 14B-only
# imaging.
# If the original mask still exists, remove it
# old_maskname = "{0}.mask".format(imagename)

# if os.path.exists(old_maskname):
#     os.system("rm -r {}".format(old_maskname))

# from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

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

# imager = PySynthesisImager(params=paramList)
imager = PyParallelContSynthesisImager(params=paramList)

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
    if restoration and niter > 0:
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

    if niter > 0:
        out_dict = imager.IBtool.getiterationsummary()

        # Save the output dictionary. Numpy should be fine for this as the
        # individual channels will get concatenated together

        np.save(imagename + ".results_dict_stage2.npy", out_dict)

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
