
'''
Create final data cubes.

Give the path to the data and which cube to create.
Puts the created cubes in a folder called "cubes" one level above where the
channel folders are.
'''

import glob
import os
import sys

from taskinit import ia

# Image product end to create cube from (i.e., image, residual, model, etc...)
output_type = sys.argv[-1]

# Name of the cube to make. The output_type is appended on
# (i.e., "M33_14B-088_HI.clean")
cube_name = sys.argv[-2]

# Path to the channel data
path = sys.argv[-3]

cube_path = os.path.join("/".join(path.rstrip("/").split("/")[:-1]), "cubes")

if not os.path.exists(cube_path):
    os.mkdir(cube_path)

search_string = "*.{}".format(output_type)
out_name = os.path.join(cube_path, "{0}.{1}".format(cube_name. output_type))

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
for chan in xrange(min_chan, max_chan + 1):
    i = 0
    while True:
        chan_num = int(channel_direcs[i].split("_")[-1])
        if chan_num == chan:
            ordered_chans.append(channel_direcs[i])
            break
        i += 1
    else:
        print("No log found for channel " + str(chan))

channel_direcs = ordered_chans

images = []
for channel in channel_direcs:
    image = glob.glob(os.path.join(channel, search_string))
    if len(image) == 0:
        print("Cannot find image in " + channel)
        continue

    images.append(image[0])

ia.imageconcat(outfile=out_name, infiles=images)
