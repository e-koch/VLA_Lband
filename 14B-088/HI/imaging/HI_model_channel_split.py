
'''
Split the HI Arecibo model for 14B-088 into individual
channels.
'''

from casa_tools import image_split_by_channel

start_chan = 0
nchan = 1231

model = "/home/ekoch/m33/14B-088/M33_14B-088_HI_model.image"
output_dir = "/home/ekoch/m33/14B-088/model_channels/"


image_split_by_channel(model, nchan=nchan, start=start_chan,
                       output_dir=output_dir)
