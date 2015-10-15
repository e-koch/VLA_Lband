
'''
Split the HI data from 14B-088 into individual channels.
'''

from casa_tools import ms_split_by_channel

vis = "/home/ekoch/m33/14B-088/14B-088_HI_LSRK.ms.contsub"
output_dir = "/home/ekoch/m33/14B-088/channel_ms/"

start_chan = 670
nchan = 210

ms_split_by_channel(vis, nchan=nchan, start=start_chan,
                    output_dir=output_dir, datacolumn='DATA')

