
'''
Split out each SPW from the combined MS (concat_and_split.py), convert
to LSRK, and subtract continuum in uv-plane
'''

import os

from tasks import mstransform

myvis = '17B-162_lines.ms'

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/17B-162/spw_setup.py"))

for spw_num in range(1, 10):

    default('mstransform')

    out_vis = "17B-162_{0}_spw_{1}_LSRK.ms.contsub"\
        .format(linespw_dict[spw_num][0], spw_num)

    casalog.post("On SPW {}".format(spw_num))

    # Note that the combined MS already only includes the calibrated data
    # with all flagged data removed.
    mstransform(vis=myvis, outputvis=out_vis, spw=spw_num, datacolumn='data',
                regridms=True, mode='channel', interpolation='fftshift',
                phasecenter='J2000 01h33m50.904 +30d39m35.79',
                restfreq=linespw_dict[spw_num][1], outframe='LSRK',
                douvcontsub=True,
                fitspw='{0}:{1}'.format(spw_num, linespw_dict[spw_num][3]),
                fitorder=0, want_cont=False)
