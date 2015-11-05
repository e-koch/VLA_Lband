
'''
Create final data cubes
'''

import glob
import os

from taskinit import ia

path = "/home/ekoch/m33/14B-088/single_channels/"
out_name = "/home/ekoch/m33/14B-088/cube_images/M33_14B-088_HI.clean.image"

channel_direcs = glob.glob(os.path.join(path, "channel_*"))

# Sort by the channel number
min_chan = 10000
max_chan = 0

for channel in channel_direcs:
    # Last 3 characters
    chan_num = int(channel.split("_")[-1])
    if chan_num < min_chan:
        min_chan = chan_num

    if chan_num > max_chan:
        max_chan = chan_num

ordered_chans = []
for chan in xrange(min_chan, max_chan+1):
    i = 0
    while True:
        chan_num = int(channel_direcs[i].split("_")[-1])
        if chan_num == chan:
            ordered_chans.append(channel_direcs[i])
            break
        i += 1
    else:
        print("No log found for channel "+str(chan))

channel_direcs = ordered_chans

images = []
for channel in channel_direcs:
    image = glob.glob(os.path.join(channel, "*.clean.image"))
    if not image:
        print("Cannot find log file in "+channel)
        continue

    images.append(image[0])

ia.imageconcat(outfile=out_name, infiles=images)
