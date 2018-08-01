

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

out_dict = tclean()

# Save the output dictionary. Numpy should be fine for this as the individual
# channels will get concatenated together

np.save(imagename + ".results_dict.npy", out_dict)
