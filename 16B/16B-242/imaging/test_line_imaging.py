
'''
Create a dirty for the given SPW

Should be run from the 16B-242_imaging folder on cedar
'''

import os
import sys
import numpy as np

# Requires analysisutils to be appended to the casa path
# Load in the auto image parameter setters
from CASA_functions import set_cellsize, set_imagesize

from tasks import tclean

spw_num = int(sys.argv[-1])

output_path = "spw_{}".format(spw_num)

if not os.path.exists(output_path):
    os.mkdir(output_path)

# Grab all of the MS tracks in the folder (should be 17)
myvis = "16B-242_lines.ms"

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/16B/spw_setup.py"))

# Assume we can set reasonable image parameters from any of the tracks
mycellsize = set_cellsize(myvis, spw_num, sample_factor=6.,
                          baseline_percentile=95,
                          return_type="str")
casalog.post("Cell size: {}".format(mycellsize))

mypblimit = 0.2

# Will look for all M33 fields and assume they are all used in the mosaic
source = 'NGC604'

myimagesize = set_imagesize(myvis, spw_num, source, sample_factor=6.,
                            max_size=15000, pblevel=mypblimit)
casalog.post("Image size: {}".format(myimagesize))

# Image ALL channels in the MS. Just looking for reduction issues
default('tclean')

# Lots of noise here, so set fairly large channels
chan_width = {0: 25, 1: 3, 2: 3, 3: 5, 4: 3,
              5: 5, 6: 5, 7: 5, 8: 3, 9: 3}

# Don't image the channel edges
pad_chan = int(np.ceil(linespw_dict[spw_num][2] * 0.05))
num_chan = int((linespw_dict[spw_num][2] - 2 * pad_chan) / chan_width[spw_num])

tclean(vis=myvis,
       datacolumn='corrected',
       imagename=os.path.join(output_path,
                              'M33_16B-242_{0}_spw_{1}.dirty'
                              .format(linespw_dict[spw_num][0], spw_num)),
       spw=str(spw_num),
       field='NGC604',
       imsize=myimagesize,
       cell=mycellsize,
       specmode='cube',
       start=pad_chan,
       width=chan_width[spw_num],
       nchan=num_chan,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       restfreq=linespw_dict[spw_num][1],
       outframe='LSRK',
       pblimit=mypblimit,
       usemask='pb',
       mask=None,
       deconvolver='hogbom',
       pbcor=False,
       veltype='radio',
       chanchunks=-1,
       restoration=False,
       parallel=True,
       )
