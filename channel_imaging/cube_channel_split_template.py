
'''
Split out channels from a cube (the mask or the model).
'''

import sys

from casa_tools import image_split_by_channel

start_chan = 0
nchan = 1231

cube = sys.argv[1]
output_dir = sys.argv[2]

image_split_by_channel(cube, nchan=nchan, start=start_chan,
                       output_dir=output_dir)
