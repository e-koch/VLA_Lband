

import sys
import os
from distutils.dir_util import mkpath

# from casa_tools import myclean
from tasks import tclean

'''
Cleans a single channel given the channel name

This is the template file. There are too many variables to easily generalize.
'''

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/17B-162/spw_setup.py"))

chan_num = int(sys.argv[-2])

wcont = sys.argv[-1]

if wcont == "T":
    myvis = "17B-162_HI_spw_0_LSRK.ms"
    out_prefix = "spw_0_perchan_wcont"
else:
    myvis = "17B-162_HI_spw_0_LSRK.ms.contsub"
    out_prefix = "spw_0_perchan"

spw_num = 0

output_path = "{0}/channel_{1}".format(out_prefix, chan_num)

# Since the imaging is not run from the parent path for the channel output,
# use `mkpath`
if not os.path.exists(output_path):
    mkpath(output_path)

mycellsize = '1.0arcsec'

myimagesize = [6250, 6144]

mypblimit = 0.1

default('tclean')

tclean(vis=myvis,
       datacolumn='corrected',
       imagename=os.path.join(output_path,
                              'M33_17B-162_{0}_channel_{1}.dirty'
                              .format(linespw_dict[spw_num][0], chan_num)),
       spw=str(spw_num),
       field='M33*',
       imsize=myimagesize,
       cell=mycellsize,
       specmode='cube',
       start=chan_num,
       width=1,
       nchan=1,
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
       parallel=False,
       )
