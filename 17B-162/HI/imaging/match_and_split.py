'''
Match and split the re-weighted HI MSs
The 14B frequency range includes 2000 channels. That same range in the 17B
data is 2006 channels. So we first regrid and split the data over the same
velocity range. The original channels are -206 and
change m/s. Regrid to something common like -210 m/s.

UPDATE: The issue with this regridding is that there is a beat pattern
when changing the channel widths by the such a small amount. CASA 5.3
*seemed* to actually be using linear interpolation, not the FFT shift.
When using the FFT shift in CASA 5.4 on other data, it caused horrific
residuals in the spectral dimension. Unsure why. However, the two
data sets have the same frequency channel width (0.977 kHz) and have a
frequency offset of 0.2% of the channel width. I'm just going to
match frequency bins to the nearest velocity instead of regridding.
For larger channel sizes, I'll round down to the nearest integer.

The individual channels are then split out for imaging. A separate
folder structure is made for each choice of channel width.

The input given is the channel width in km/s. It is assumed that
the start and end points will be the same (or rounded up by 1),
'''

import numpy as np
import sys
from glob import glob
import os

from tasks import mstransform

# Use astropy's spectral conversion
# Needs to be installed separately
from astropy import units as u

use_contsub = True if sys.argv[-2] == "True" else False

# All in km/s
chan_width = float(sys.argv[-1])
start_vel = -330
end_vel = -50

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

# Old version

# seventeenB_ms_regrid = "{0}.{1}.regrid".format(seventeenB_ms, chan_width_label)

# if os.path.exists(seventeenB_ms_regrid):
#     casalog.post("Found the regridded 17B MS. Skipping mstransform.")
# else:
#     casalog.post("Regridding 17B")

#     mstransform(vis=seventeenB_ms,
#                 outputvis=seventeenB_ms_regrid,
#                 datacolumn='data',
#                 spw='0', regridms=True, mode='velocity', veltype='radio',
#                 start="{}km/s".format(start_vel),
#                 width="{}km/s".format(chan_width), nchan=n_chan,
#                 interpolation='fftshift', restfreq='1.420405752GHz',
#                 createmms=True, separationaxis='auto', numsubms=31)

# # Now the 14B data. Only keep the fields used in the 17B data
# fourteenB_ms_regrid = "{0}.{1}.regrid".format(fourteenB_ms, chan_width_label)

# if os.path.exists(fourteenB_ms_regrid):
#     casalog.post("Found the regridded 14B MS. Skipping mstransform.")
# else:

#     casalog.post("Regridding 14B")

#     # Also convert th 14B data to an MMS w/ the same number of sub-MS as 17B
#     mstransform(vis=fourteenB_ms,
#                 outputvis=fourteenB_ms_regrid,
#                 datacolumn='data',
#                 field='M33_2,M33_6,M33_7_center,M33_8,M33_11,M33_12,M33_14',
#                 spw='0', regridms=True, mode='velocity', veltype='radio',
#                 start="{}km/s".format(start_vel),
#                 width="{}km/s".format(chan_width), nchan=n_chan,
#                 interpolation='fftshift', restfreq='1.420405752GHz',
#                 createmms=True, separationaxis='auto', numsubms=31)

# New version:


def vel_to_freq(vel_or_freq, rest_freq=1.42040575177 * u.GHz,
                unit=u.Hz):
    '''
    Using radio velocity here.
    '''
    equiv = u.doppler_radio(rest_freq)

    return vel_or_freq.to(unit, equiv)


def closest_channel(freqs, targ_freq):
    return np.argmin(np.abs(freqs - targ_freq))


#n_chan = int(np.ceil((end_vel - start_vel) / chan_width))

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

if start == nchan:
    casalog.post("No more channels to split!")
    raise ValueError("No more channels to split!")

for chan in range(start, nchan + 1):

    casalog.post("On splitting channel {}".format(chan))

    ind_chan_path = os.path.join(chan_path,
                                 "channel_{}".format(chan))
    if not os.path.exists(ind_chan_path):
        os.mkdir(ind_chan_path)

    sevenB_split_msname = "{0}_channel_{1}".format(seventeenB_ms, chan)

    start_17B = chan * navg_channel + start_17B_chan
    end_17B = (chan + 1) * navg_channel + start_17B_chan - 1

    mstransform(vis=seventeenB_ms,
                outputvis=os.path.join(ind_chan_path, sevenB_split_msname + ".ms"),
                datacolumn='data',
                mode='channel',
                field=myfields,
                spw='0:{0}~{1}'.format(start_17B, end_17B),
                chanaverage=True,
                chanbin=navg_channel)

    fourB_split_msname = "{0}_channel_{1}.ms".format(fourteenB_ms, chan)

    start_14B = chan * navg_channel + start_14B_chan
    end_14B = (chan + 1) * navg_channel + start_14B_chan - 1

    if navg_channel == 1:
        # These should equal when not averaging channels
        assert start_14B == end_14B
        spw_selec = "0:{0}".format(start_14B)
    else:
        spw_selec = '0:{0}~{1}'.format(start_14B, end_14B)

    mstransform(vis=fourteenB_ms,
                outputvis=os.path.join(ind_chan_path, fourB_split_msname),
                datacolumn='data',
                mode='channel',
                field=myfields,
                spw=spw_selec,
                chanaverage=True if navg_channel > 1 else False,
                chanbin=navg_channel)
