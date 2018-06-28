
'''
Create a dirty for the given SPW

Should be run from the 17B-162_imaging folder on cedar
'''

import os
import sys
from glob import glob

# Requires analysisutils to be appended to the casa path
# Load in the auto image parameter setters
from CASA_functions import set_cellsize, set_imagesizes

from tasks import tclean, impbcor

spw_num = sys.argv[-1]

output_path = "spw_{}".format(spw_num)

if not os.path.exists(output_path):
    os.mkdir(output_path)

# Grab all of the MS tracks in the folder (should be 17)
myvis = glob("*.speclines.ms")

assert len(myvis) == 17

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/17B-162/spw_setup.py"))

# Assume we can set reasonable image parameters from any of the tracks
mycellsize = set_cellsize(myvis[0], spw_num, sample_factor=6.,
                          baseline_percentile=95,
                          return_type="str")
print("Cell size: {}".format(mycellsize))

mypblimit = 0.1

# Will look for all M33 fields and assume they are all used in the mosaic
source = 'M33'

myimagesize = set_imagesize(myvis[0], spw_num, source, sample_factor=6.,
                            pblevel=0.1, max_size=15000, pblevel=mypblimit)

# Image ALL channels in the MS. Just looking for reduction issues
default('tclean')

tclean(vis=myvis,
       datacolumn='corrected',
       imagename=os.path.join(output_path,
                              'M33_14B-088_{0}_spw_{1}.dirty'.format(linespw_dict[spw_num][0],
                                                                     spw_num)),
       field='M33*',
       imsize=myimagesize,
       cell=mycellsize,
       specmode='cube',
       start=1,
       width=1,
       nchan=-1,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h33m50.904 +30d39m35.79',
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
