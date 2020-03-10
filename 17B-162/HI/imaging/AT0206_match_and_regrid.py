
'''
Try including the old AT0206 data with the 14B and 17B. The data is in the
old VLA format with individual SPWs for each correlation. It needs to be
regridded.

The 14B and 17B data will have to match the AT0206 velocity extents, as well.
The velocity extents are -64.2 to -299.9 km/s. The velocity channels have a
width of 1.288 km/s. May as well regrid to 1.3 km/s.

Potential issues:

* Regridding into a single SPW

* Different pointing centres

'''

# from paths import
import os
from glob import glob
import numpy as np

from tasks import mstransform

mstransform(vis='M33_AT0206_b_c_LSRK.ms',
            outputvis="M33_AT0206_b_c_LSRK_combSPW.ms",
            regridms=True, combinespws=True, datacolumn='data',
            phasecenter="J2000 01h33m50.904 +30d39m35.79",
            restfreq="1420.40575177MHz", mode='velocity',
            nchan=-1, start='-299km/s', width='1.3km/s',
            interpolation='fftshift', outframe='LSRK',
            )

chan_width = 1.3

n_chan = 192

chan_width_label = "{}kms".format(chan_width).replace(".", "")

chan_path = "HI_{0}_{1}".format("contsub",
                                chan_width_label)

if not os.path.exists(chan_path):

    os.mkdir(chan_path)

    start = 0

else:

    # Check if some of the channels have already been split
    exist_chans = glob("{}/channel_*".format(chan_path))

    # Didn't find any existing channels
    if len(exist_chans) == 0:
        start = 0

    else:
        # Pick out the channel numbers
        nums = np.array([int(chan.split("_")[-1]) for chan in exist_chans])

        # Assume that the last job stopped while part of the way through
        # the final channel
        final_chan = exist_chans[nums.argmax()]

        os.system("rm -r {}".format(final_chan))

        start = nums.max() - 1

if start == n_chan:
    casalog.post("No more channels to split!")
    raise ValueError("No more channels to split!")

for chan in range(start, n_chan + 1):

    casalog.post("On splitting channel {}".format(chan))

    ind_chan_path = os.path.join(chan_path,
                                 "channel_{}".format(chan))

    split_msname = "{0}_channel_{1}".format(seventeenB_ms_regrid, chan)

    if not os.path.exists(ind_chan_path):
        os.mkdir(ind_chan_path)

    mstransform(vis=seventeenB_ms_regrid,
                outputvis=os.path.join(ind_chan_path, split_msname + ".mms"),
                datacolumn='data',
                field='M33*',
                spw='0:{}'.format(chan))

    # You know what has a lot of files? MMSs! The cluster doesn't
    # like this, so convert back to an MS
    default('split')
    split(vis=os.path.join(ind_chan_path, split_msname + ".mms"),
          outputvis=os.path.join(ind_chan_path, split_msname), keepmms=False,
          datacolumn='data')

    os.system("rm -r {}".format(os.path.join(ind_chan_path, split_msname + ".mms")))

    split_14B_msname = "{0}_channel_{1}".format(fourteenB_ms_regrid,
                                                chan)

    mstransform(vis=fourteenB_ms_regrid,
                outputvis=os.path.join(ind_chan_path, split_14B_msname + ".mms"),
                datacolumn='data',
                field='M33*',
                spw='0:{}'.format(chan))

    split(vis=os.path.join(ind_chan_path, split_14B_msname + ".mms"),
          outputvis=os.path.join(ind_chan_path, split_14B_msname), keepmms=False,
          datacolumn='data')

    os.system("rm -r {}".format(os.path.join(ind_chan_path, split_14B_msname + ".mms")))
