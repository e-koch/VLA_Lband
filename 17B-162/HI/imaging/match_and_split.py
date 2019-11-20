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

from tasks import mstransform, partition, split, concat

# Use astropy's spectral conversion
# Needs to be installed separately
from astropy import units as u

# This is here for local runs to avoid needing to make an MMS
# Mostly for storage reasons.
use_parallel = True

use_contsub = True if sys.argv[-5] == "True" else False

# All in km/s
chan_width = float(sys.argv[-4])
start_vel = -330
end_vel = -50

part = int(sys.argv[-3])
total_parts = int(sys.argv[-2])

scratch_path = str(sys.argv[-1])

chan_width_label = "{}kms".format(chan_width).replace(".", "_")
chan_width_quant = chan_width * u.km / u.s

chan_path = "HI_{0}_{1}".format("contsub" if use_contsub else "nocontsub",
                                chan_width_label)

out_path = os.path.join(scratch_path, chan_path)

# Common fields in B and C
myfields = 'M33_2,M33_6,M33_7_center,M33_8,M33_11,M33_12,M33_14'

# ~1334 for 0.21 km/s
# n_chan = int(np.ceil((end_vel - start_vel) / chan_width))


chunk_suffix = "{0}.mms_chunk_{1}".format(chan_width_label, part)


if use_contsub:

    fourteenB_ms_orig = "{0}/14B-088_HI_LSRK.ms.contsub".format(scratch_path)
    seventeenB_ms_orig = "{0}/17B-162_HI_spw_0_LSRK.ms.contsub".format(scratch_path)

    fourteenB_ms = "14B-088_HI_LSRK.ms.contsub.{0}".format(chunk_suffix)
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.ms.contsub.{0}".format(chunk_suffix)

    concat_vis = "14B_17B_HI_LSRK.ms.contsub"

else:

    fourteenB_ms_orig = "{0}/14B-088_HI_LSRK.ms".format(scratch_path)
    seventeenB_ms_orig = "{0}/17B-162_HI_spw_0_LSRK.ms".format(scratch_path)

    fourteenB_ms = "14B-088_HI_LSRK.ms.{0}".format(chunk_suffix)
    seventeenB_ms = "17B-162_HI_spw_0_LSRK.ms.{0}".format(chunk_suffix)

    concat_vis = "14B_17B_HI_LSRK.ms"


# The chunks are already an MMS. So don't worry about this part.
seventeenB_mms = seventeenB_ms
fourteenB_mms = fourteenB_ms

# Create an MMS prior to splitting to that the split can be run in parallel

# if use_parallel:

#     if os.path.exists(seventeenB_mms):
#         casalog.post("Found the 17B MMS. Skipping mstransform.")
#     else:
#         casalog.post("Regridding 17B")

#         partition(vis=seventeenB_ms,
#                   outputvis=seventeenB_mms,
#                   createmms=True,
#                   flagbackup=False,
#                   numsubms=31)

#         # mstransform(vis=seventeenB_ms,
#         #             outputvis=seventeenB_mms,
#         #             datacolumn='data',
#         #             spw='0', regridms=True, mode='velocity', veltype='radio',
#         #             start="{}km/s".format(start_vel),
#         #             width="{}km/s".format(chan_width), nchan=n_chan,
#         #             interpolation='fftshift', restfreq='1.420405752GHz',
#         #             createmms=True, separationaxis='auto', numsubms=31)

#     # Now the 14B data. Only keep the fields used in the 17B data

#     if os.path.exists(fourteenB_mms):
#         casalog.post("Found the regridded 14B MS. Skipping mstransform.")
#     else:

#         casalog.post("Regridding 14B")

#         partition(vis=fourteenB_ms,
#                   outputvis=fourteenB_mms,
#                   createmms=True,
#                   flagbackup=False,
#                   numsubms=31)

#     #     # Also convert th 14B data to an MMS w/ the same number of sub-MS as 17B
#     #     mstransform(vis=fourteenB_ms,
#     #                 outputvis=fourteenB_ms_regrid,
#     #                 datacolumn='data',
#     #                 field='M33_2,M33_6,M33_7_center,M33_8,M33_11,M33_12,M33_14',
#     #                 spw='0', regridms=True, mode='velocity', veltype='radio',
#     #                 start="{}km/s".format(start_vel),
#     #                 width="{}km/s".format(chan_width), nchan=n_chan,
#     #                 interpolation='fftshift', restfreq='1.420405752GHz',
#     #                 createmms=True, separationaxis='auto', numsubms=31)
# else:

#     seventeenB_mms = seventeenB_ms
#     fourteenB_mms = fourteenB_ms

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
tb.open(os.path.join(fourteenB_ms_orig, 'SPECTRAL_WINDOW'))
chanfreqs_14B = tb.getcol('CHAN_FREQ').squeeze()
tb.close()

delta_freq_14B = np.abs(np.diff(chanfreqs_14B))[0]


tb.open(os.path.join(seventeenB_ms_orig, 'SPECTRAL_WINDOW'))
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

if not os.path.exists(chan_path):

    os.mkdir(chan_path)

    start = 0

# else:

#     # Check if some of the channels have already been split
#     exist_chans = glob("{}/channel_*".format(chan_path))

#     # Didn't find any existing channels
#     if len(exist_chans) == 0:
#         start = 0

#     else:
#         # Pick out the channel numbers
#         nums = np.array([int(chan.split("_")[-1]) for chan in exist_chans])

#         # Assume that the last job stopped while part of the way through
#         # the final channel
#         final_chan = exist_chans[nums.argmax()]

#         os.system("rm -r {}".format(final_chan))

#         start = nums.max() - 1

# if start == nchan:
#     casalog.post("No more channels to split!")
#     raise ValueError("No more channels to split!")

nchan_part = int(np.ceil(nchan / total_parts))

start = part * nchan_part
end = min((part + 1) * nchan_part, nchan)

# for chan in range(start, nchan + 1):
for chan_chunk, chan in enumerate(range(start, end)):

    casalog.post("On splitting channel {}".format(chan))

    ind_chan_path = os.path.join(chan_path,
                                 "channel_{}".format(chan))
    if not os.path.exists(ind_chan_path):
        os.mkdir(ind_chan_path)

    # Does the concat MS exist already? If yes, skip
    concat_vis_name = '14B_17B_channel_{}.ms'.format(chan)

    concat_ms = os.path.join(ind_chan_path, concat_vis_name)
    if os.path.exists(concat_ms):
        casalog.post("Channel {} already has an concatenated MS. Skipping.".format(chan))
        continue

    if use_parallel:
        out_channel = os.path.join(out_path, "channel_{}".format(chan))
        if not os.path.exists(out_channel):
            os.mkdir(out_channel)
        else:
            # Does the MS already exist there? If so, skip it here.
            scratch_ms = os.path.join(out_channel, concat_vis_name)
            if os.path.exists(scratch_ms):
                casalog.post("Found the split + concat MS for {} in scratch. "
                             "Skipping.".format(chan))
                continue

    sevenB_split_msname = "{0}_channel_{1}.ms".format(seventeenB_ms, chan)
    if use_parallel:
        sevenB_split_mmsname = "{0}_channel_{1}.mms".format(seventeenB_ms,
                                                            chan)
    else:
        sevenB_split_mmsname = sevenB_split_msname

    # Because of the chunks, we now just iterate through without the offsets
    start_17B = chan_chunk * navg_channel
    end_17B = (chan_chunk + 1) * navg_channel + - 1

    # start_17B = chan * navg_channel + start_17B_chan
    # end_17B = (chan + 1) * navg_channel + start_17B_chan - 1

    if navg_channel == 1:
        # These should equal when not averaging channels
        assert start_17B == end_17B
        spw_selec = "0:{0}".format(start_17B)
    else:
        spw_selec = '0:{0}~{1}'.format(start_17B, end_17B)

    mstransform(vis=seventeenB_mms,
                outputvis=os.path.join(ind_chan_path, sevenB_split_mmsname),
                datacolumn='data',
                mode='channel',
                field=myfields,
                spw=spw_selec,
                chanaverage=True if navg_channel > 1 else False,
                chanbin=navg_channel)

    fourB_split_msname = "{0}_channel_{1}.ms".format(fourteenB_ms, chan)
    if use_parallel:
        fourB_split_mmsname = "{0}_channel_{1}.mms".format(fourteenB_ms, chan)
    else:
        fourB_split_mmsname = fourB_split_msname

    start_14B = chan_chunk * navg_channel
    end_14B = (chan_chunk + 1) * navg_channel + - 1

    # start_14B = chan * navg_channel + start_14B_chan
    # end_14B = (chan + 1) * navg_channel + start_14B_chan - 1

    if navg_channel == 1:
        # These should equal when not averaging channels
        assert start_14B == end_14B
        spw_selec = "0:{0}".format(start_14B)
    else:
        spw_selec = '0:{0}~{1}'.format(start_14B, end_14B)

    mstransform(vis=fourteenB_mms,
                outputvis=os.path.join(ind_chan_path, fourB_split_mmsname),
                datacolumn='data',
                mode='channel',
                field=myfields,
                spw=spw_selec,
                chanaverage=True if navg_channel > 1 else False,
                chanbin=navg_channel)

    if use_parallel:
        # Convert the final MMS to an MS b/c an MMS uses a lot of files and
        # clusters don't like that.
        split(vis=os.path.join(ind_chan_path, sevenB_split_mmsname),
              outputvis=os.path.join(ind_chan_path, sevenB_split_msname),
              keepmms=False, datacolumn='DATA')

        split(vis=os.path.join(ind_chan_path, fourB_split_mmsname),
              outputvis=os.path.join(ind_chan_path, fourB_split_msname),
              keepmms=False, datacolumn='DATA')

        # Remove per-chan MMS
        os.system("rm -rf {}".format(os.path.join(ind_chan_path,
                                                  sevenB_split_mmsname)))
        os.system("rm -rf {}".format(os.path.join(ind_chan_path,
                                                  fourB_split_mmsname)))

    # tclean calls have been ignoring the C-config data
    # concat the channel MSs for imaging

    # If the concat ms already exists, delete it. Otherwise more data
    # will be appended on
    if os.path.exists(concat_ms):
        os.system("rm -rf {}".format(concat_ms))

    concat(vis=[os.path.join(ind_chan_path, sevenB_split_msname),
                os.path.join(ind_chan_path, fourB_split_msname)],
           concatvis=concat_ms)

    # Clean-up temporary MS
    os.system("rm -rf {}".format(os.path.join(ind_chan_path,
                                              sevenB_split_msname)))
    os.system("rm -rf {}".format(os.path.join(ind_chan_path,
                                              fourB_split_msname)))
