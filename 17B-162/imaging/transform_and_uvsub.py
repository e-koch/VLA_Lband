
'''
Split out each SPW from the combined MS (concat_and_split.py), convert
to LSRK, and subtract continuum in uv-plane
'''

import os

from tasks import mstransform, uvcontsub, partition, split

myvis = '17B-162_lines.ms'

# Load in the SPW dict in the repo on cedar
execfile(os.path.expanduser("~/code/VLA_Lband/17B-162/spw_setup.py"))

for spw_num in range(0, 10):

    default('mstransform')

    casalog.post("On SPW {}".format(spw_num))

    # Note that the combined MS already only includes the calibrated data
    # with all flagged data removed.

    # If this is HI, we want to keep the continuum version to look for
    # absorption features. Split the uvsubtraction into a separate function

    if spw_num == 0:
        do_uvcontsub = False
        out_vis = "17B-162_{0}_spw_{1}_LSRK.ms"\
            .format(linespw_dict[spw_num][0], spw_num)
    else:
        do_uvcontsub = True
        out_vis = "17B-162_{0}_spw_{1}_LSRK.ms.contsub"\
            .format(linespw_dict[spw_num][0], spw_num)

    mstransform(vis=myvis, outputvis=out_vis, spw=str(spw_num), datacolumn='data',
                regridms=True, mode='channel', interpolation='fftshift',
                phasecenter='J2000 01h33m50.904 +30d39m35.79',
                restfreq=linespw_dict[spw_num][1], outframe='LSRK',
                douvcontsub=do_uvcontsub,
                fitspw='{0}:{1}'.format(spw_num, linespw_dict[spw_num][3]),
                fitorder=0, want_cont=False)

    # Separate uvcontsub for HI
    if spw_num == 0:
        out_vis_cs = "17B-162_{0}_spw_{1}_LSRK.ms.contsub"\
            .format(linespw_dict[spw_num][0], spw_num)

        out_mms_vis_cs = "17B-162_{0}_spw_{1}_LSRK.mms.contsub"\
             .format(linespw_dict[spw_num][0], spw_num)

        # The operation is much fast in parallel, so make an MMS and then
        # convert back
        partition(vis=out_vis, outputvis=out_vis[:-3] + ".mms", createmms=True,
                  separationaxis='auto', flagbackup=False)

        uvcontsub(vis=out_vis[:-3] + ".mms",
                  fitspw='{0}:{1}'.format(spw_num, linespw_dict[spw_num][3]),
                  fitorder=0, want_cont=False)

        default('split')
        split(vis=out_mms_vis_cs, outputvis=out_vis_cs, keepmms=False)
