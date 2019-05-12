
'''
Split out each SPW from the combined MS (concat_and_split.py), convert
to LSRK, and subtract continuum in uv-plane
'''

import os
import sys

from tasks import mstransform, uvcontsub

myvis = '17B-162_lines.ms'

spw_num = int(sys.argv[-1])

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/17B-162/spw_setup.py"))


default('mstransform')

casalog.post("On SPW {}".format(spw_num))

# Note that the combined MS already only includes the calibrated data
# with all flagged data removed.

out_vis = "17B-162_{0}_spw_{1}_LSRK.ms"\
    .format(linespw_dict[spw_num][0], spw_num)


# in casa 5.4.1, fftshift does something that is not linear interpolation
# but it gives severe edge effects when I've used it on the 13A-213 data
mstransform(vis=myvis, outputvis=out_vis, spw=str(spw_num),
            datacolumn='data',
            field='M33*',
            regridms=True, mode='channel',
            interpolation='linear',  # 'fftshift',
            phasecenter='J2000 01h33m50.904 +30d39m35.79',
            restfreq=linespw_dict[spw_num][1], outframe='LSRK',
            douvcontsub=False)

default('uvcontsub')

# Separate uvcontsub call
uvcontsub(vis=out_vis,
          fitspw='{0}:{1}'.format(spw_num, linespw_dict[spw_num][3]),
          fitorder=0, want_cont=False)