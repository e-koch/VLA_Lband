
'''
Match and combine the re-weighted HI MSs

The 14B frequency range includes 2000 channels. That same range in the 17B
data is 2006 channels. So we first regrid and split the data over the same
velocity range before concatenating. The original channels are -206 and
change m/s. Regrid to something common like -210 m/s
'''

import numpy as np

from tasks import mstransform, virtualconcat, rmtables

use_contsub = True

if use_contsub:
    fourteenB_ms = "14B-088_HI_LSRK.ms.contsub"
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.mms.contsub"

    concat_vis = "14B_17B_HI_LSRK.mms.contsub"

else:

    fourteenB_ms = "14B-088_HI_LSRK.ms"
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.mms"

    concat_vis = "14B_17B_HI_LSRK.mms"

# All in km/s
chan_width = 0.21
start_vel = -330
end_vel = -50

# ~1334
n_chan = int(np.ceil((end_vel - start_vel) / chan_width))

seventeenB_ms_regrid = "{}.regrid".format(seventeenB_ms)

mstransform(vis=seventeenB_ms,
            outputvis=seventeenB_ms_regrid,
            spw='0', regridms=True, mode='velocity', veltype='radio',
            start="{}km/s".format(start_vel),
            width="{}km/s".format(chan_width), nchan=n_chan,
            interpolation='fftshift', restfreq='1.420405752GHz')

# Now the 14B data. Only keep the fields used in the 17B data
fourteenB_ms_regrid = "{}.regrid".format(fourteenB_ms)

# Also convert th 14B data to an MMS w/ the same number of sub-MS as 17B
mstransform(vis=fourteenB_ms,
            outputvis=fourteenB_ms_regrid,
            field='M33_2,M33_6,M33_7_center,M33_8,M33_11,M33_12,M33_14',
            spw='0', regridms=True, mode='velocity', veltype='radio',
            start="{}km/s".format(start_vel),
            width="{}km/s".format(chan_width), nchan=n_chan,
            interpolation='fftshift', restfreq='1.420405752GHz',
            createmms=True, separationaxis='auto', numsubms=31)

# concatenate together, and will also remove the individual MS
virtualconcat(vis=[fourteenB_ms_regrid, seventeenB_ms_regrid],
              concatvis=concat_vis)
