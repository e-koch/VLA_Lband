
'''
The 17B MS is 1.4 TB and that's too big to move to the SSD storage on the
cedar nodes. Instead, I'm going to make chunks of the 14B and 17B MSs that
`match_and_split.py` will use.
'''


import numpy as np
import sys
from glob import glob
import os

from tasks import mstransform

# Use astropy's spectral conversion
# Needs to be installed separately
from astropy import units as u

# This is here for local runs to avoid needing to make an MMS
# Mostly for storage reasons.
use_parallel = True

use_contsub = True if sys.argv[-3] == "True" else False

# All in km/s
chan_width = float(sys.argv[-2])
start_vel = -330
end_vel = -50

total_parts = int(sys.argv[-1])

chan_width_label = "{}kms".format(chan_width).replace(".", "_")
chan_width_quant = chan_width * u.km / u.s

# Common fields in B and C
myfields = 'M33_2,M33_6,M33_7_center,M33_8,M33_11,M33_12,M33_14'

# ~1334 for 0.21 km/s
# n_chan = int(np.ceil((end_vel - start_vel) / chan_width))

if use_contsub:
    fourteenB_ms = "14B-088_HI_LSRK.ms.contsub"
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.ms.contsub"

    concat_vis = "14B_17B_HI_LSRK.ms.contsub"

else:

    fourteenB_ms = "14B-088_HI_LSRK.ms"
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.ms"

    concat_vis = "14B_17B_HI_LSRK.ms"


def vel_to_freq(vel_or_freq, rest_freq=1.42040575177 * u.GHz,
                unit=u.Hz):
    '''
    Using radio velocity here.
    '''
    equiv = u.doppler_radio(rest_freq)

    return vel_or_freq.to(unit, equiv)


def closest_channel(freqs, targ_freq):
    return np.argmin(np.abs(freqs - targ_freq))



# Get the HI SPW freqs
tb.open(os.path.join(fourteenB_ms, 'SPECTRAL_WINDOW'))
chanfreqs_14B = tb.getcol('CHAN_FREQ').squeeze()
tb.close()

delta_freq_14B = np.abs(np.diff(chanfreqs_14B))[0]


tb.open(os.path.join(seventeenB_ms, 'SPECTRAL_WINDOW'))
chanfreqs_17B = tb.getcol('CHAN_FREQ').squeeze()
tb.close()

delta_freq_17B = np.abs(np.diff(chanfreqs_17B))[0]

# They should be really close
assert abs(delta_freq_14B - delta_freq_17B) < 0.1

# Find the number of channels to get closest to the requested velocity width
vunit = u.km / u.s
vel_width = \
    np.abs(vel_to_freq(chanfreqs_14B[len(chanfreqs_14B) // 2] * u.Hz, unit=vunit) -
           vel_to_freq(chanfreqs_14B[len(chanfreqs_14B) // 2 - 1] * u.Hz, unit=vunit))

navg_channel = int(round((chan_width_quant / vel_width).value))

start_freq = vel_to_freq(start_vel * u.km / u.s)
end_freq = vel_to_freq(end_vel * u.km / u.s)

# Find the start and finish channels in each MS
start_14B_chan = closest_channel(chanfreqs_14B * u.Hz, start_freq)
end_14B_chan = closest_channel(chanfreqs_14B * u.Hz, end_freq)

if start_14B_chan > end_14B_chan:
    start_14B_chan, end_14B_chan = end_14B_chan, start_14B_chan

start_17B_chan = closest_channel(chanfreqs_17B * u.Hz, start_freq)
end_17B_chan = closest_channel(chanfreqs_17B * u.Hz, end_freq)

if start_17B_chan > end_17B_chan:
    start_17B_chan, end_17B_chan = end_17B_chan, start_17B_chan

# Channel number in terms of original channel widths
nchan_14B = end_14B_chan - start_14B_chan
nchan_17B = end_17B_chan - start_17B_chan

# These need to be the same. Catch possible rounding errors
assert nchan_14B == nchan_17B

# Now convert to the number of channels at the expected velocity resolution
nchan = nchan_14B // navg_channel
# Pad number to reach integer factor of navg_channel
if nchan_14B % navg_channel != 0:
    nchan += 1

nchan_part = int(np.ceil(nchan / total_parts))


for part in range(total_parts):

    start = part * nchan_part
    end = min((part + 1) * nchan_part, nchan)

    casalog.post("On splitting chunk {}".format(part))

    seventeenB_mms = "{0}.{1}.mms_chunk_{2}".format(seventeenB_ms,
                                                    chan_width_label, part)
    fourteenB_mms = "{0}.{1}.mms_chunk_{2}".format(fourteenB_ms,
                                                   chan_width_label, part)

    start_17B = start + start_17B_chan
    end_17B = (end + 1) + start_17B_chan - 1

    spw_selec_17B = "0:{0}~{1}".format(start_17B, end_17B)

    mstransform(vis=seventeenB_ms,
                outputvis=seventeenB_mms,
                datacolumn='data',
                mode='channel',
                field=myfields,
                spw=spw_selec_17B,
                chanaverage=False,
                createmms=True,
                separationaxis='auto',
                numsubms=31)

    start_14B = start + start_14B_chan
    end_14B = (end + 1) + start_14B_chan - 1

    spw_selec_14B = "0:{0}~{1}".format(start_14B, end_14B)

    mstransform(vis=fourteenB_ms,
                outputvis=fourteenB_mms,
                datacolumn='data',
                mode='channel',
                field=myfields,
                spw=spw_selec_14B,
                chanaverage=False,
                createmms=True,
                separationaxis='auto',
                numsubms=31)
