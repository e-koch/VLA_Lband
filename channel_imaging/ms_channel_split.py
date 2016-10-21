
'''
Split the HI data from 14B-088 into individual channels.
'''

import sys
import os

from casa_tools import ms_split_by_channel

vis = sys.argv[-6]
output_dir = sys.argv[-5]
spw = sys.argv[-4]
field = sys.argv[-3]

# Create folder for the split channel MSs
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# start_chan = 670
# nchan = 210
start_chan = int(sys.argv[-2])
nchan = int(sys.argv[-1])

ms_split_by_channel(vis, nchan=nchan, start=start_chan,
                    output_dir=output_dir, datacolumn='DATA',
                    field=field, spw=spw, use_split=True)
