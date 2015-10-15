
'''
Split the HI Arecibo mask for 14B-088 into individual
channels.
'''

from casa_tools import image_split_by_channel

start_chan = 670
nchan = 400

mask = "/home/ekoch/m33/14B-088/M33_14B-088_HI_mask.image"
output_dir = "/home/ekoch/m33/14B-088/mask_channels/"


image_split_by_channel(mask, nchan=nchan, start=start_chan,
                       output_dir=output_dir)
