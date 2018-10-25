
'''
Split out each SPW from the combined MS (concat_and_split.py), convert
to LSRK, and subtract continuum in uv-plane
'''

import os
import sys

from tasks import mstransform, uvcontsub, partition, split

myvis = '16B-236_lines.ms'

spw_num = int(sys.argv[-1])

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/16B/spw_setup.py"))


default('mstransform')

casalog.post("On SPW {}".format(spw_num))

# Note that the combined MS already only includes the calibrated data
# with all flagged data removed.

# If this is HI, we want to keep the continuum version to look for
# absorption features. Split the uvsubtraction into a separate function

out_vis = "16B-236_{0}_spw_{1}_LSRK.ms"\
    .format(linespw_dict[spw_num][0], spw_num)

out_vis_mms = "16B-236_{0}_spw_{1}_LSRK.mms"\
    .format(linespw_dict[spw_num][0], spw_num)

mstransform(vis=myvis, outputvis=out_vis_mms, spw=str(spw_num),
            datacolumn='data',
            regridms=True, mode='channel', interpolation='fftshift',
            # phasecenter='J2000 01h33m50.904 +30d39m35.79',
            restfreq=linespw_dict[spw_num][1], outframe='LSRK',
            douvcontsub=False)

# Separate uvcontsub

out_vis_cs = "16B-236_{0}_spw_{1}_LSRK.ms.contsub"\
    .format(linespw_dict[spw_num][0], spw_num)

out_mms_vis_cs = "16B-236_{0}_spw_{1}_LSRK.mms.contsub"\
     .format(linespw_dict[spw_num][0], spw_num)

# The operation is much fast in parallel, so make an MMS and then
# convert back
# default('partition')
# partition(vis=out_vis, outputvis=out_vis[:-3] + ".mms", createmms=True,
#           separationaxis='auto', flagbackup=False)

default('uvcontsub')
uvcontsub(vis=out_vis_mms,
          fitspw='0:{1}'.format(linespw_236_uvsub[spw_num]),
          fitorder=0 if spw_num != 0 else 1,
          want_cont=False)

default('split')
split(vis=out_mms_vis_cs, outputvis=out_vis_cs, keepmms=False,
      datacolumn='data')
os.system("rm -r {}".format(out_mms_vis_cs))

default('split')
split(vis=out_vis_mms, outputvis=out_vis, keepmms=False,
      datacolumn='data')

os.system("rm -r {}".format(out_vis_mms))
