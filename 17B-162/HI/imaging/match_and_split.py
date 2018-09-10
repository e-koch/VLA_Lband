'''
Match and split the re-weighted HI MSs
The 14B frequency range includes 2000 channels. That same range in the 17B
data is 2006 channels. So we first regrid and split the data over the same
velocity range. The original channels are -206 and
change m/s. Regrid to something common like -210 m/s.

The individual channels are then split out for imaging. A separate
folder structure is made for each choice of channel width.

The input given is the channel width in km/s. It is assumed that
the start and end points will be the same (or rounded up by 1),
'''

import numpy as np
import sys
from glob import glob
import os

from tasks import mstransform, virtualconcat, rmtables

use_contsub = True if sys.argv[-2] == "True" else False

# All in km/s
chan_width = float(sys.argv[-1])
start_vel = -330
end_vel = -50

chan_width_label = "{}kms".format(chan_width).replace(".", "")

# ~1334 for 0.21 km/s
n_chan = int(np.ceil((end_vel - start_vel) / chan_width))

if use_contsub:
    fourteenB_ms = "14B-088_HI_LSRK.ms.contsub"
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.ms.contsub"

    concat_vis = "14B_17B_HI_LSRK.ms.contsub"

else:

    fourteenB_ms = "14B-088_HI_LSRK.ms"
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.ms"

    concat_vis = "14B_17B_HI_LSRK.ms"

seventeenB_ms_regrid = "{0}.{1}.regrid".format(seventeenB_ms, chan_width_label)

if os.path.exists(seventeenB_ms_regrid):
    casalog.post("Found the regridded 17B MS. Skipping mstransform.")
else:
    casalog.post("Regridding 17B")

    mstransform(vis=seventeenB_ms,
                outputvis=seventeenB_ms_regrid,
                datacolumn='data',
                spw='0', regridms=True, mode='velocity', veltype='radio',
                start="{}km/s".format(start_vel),
                width="{}km/s".format(chan_width), nchan=n_chan,
                interpolation='fftshift', restfreq='1.420405752GHz',
                createmms=True, separationaxis='auto', numsubms=31)

# Now the 14B data. Only keep the fields used in the 17B data
fourteenB_ms_regrid = "{0}.{1}.regrid".format(fourteenB_ms, chan_width_label)

if os.path.exists(fourteenB_ms_regrid):
    casalog.post("Found the regridded 14B MS. Skipping mstransform.")
else:

    casalog.post("Regridding 14B")

    # Also convert th 14B data to an MMS w/ the same number of sub-MS as 17B
    mstransform(vis=fourteenB_ms,
                outputvis=fourteenB_ms_regrid,
                datacolumn='data',
                field='M33_2,M33_6,M33_7_center,M33_8,M33_11,M33_12,M33_14',
                spw='0', regridms=True, mode='velocity', veltype='radio',
                start="{}km/s".format(start_vel),
                width="{}km/s".format(chan_width), nchan=n_chan,
                interpolation='fftshift', restfreq='1.420405752GHz',
                createmms=True, separationaxis='auto', numsubms=31)

# Now split out individual channels for imaging.

chan_path = "HI_{0}_{1}".format("contsub" if use_contsub else "nocontsub",
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
