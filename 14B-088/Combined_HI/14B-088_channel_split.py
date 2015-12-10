
'''
Split the regridded HI data from 14B-088 into individual channels.
'''

from casa_tools import ms_split_by_channel

vis = "/global/scratch/ekoch/combined/14B-088_HI_LSRK_AT0206_regrid.ms.contsub"
output_dir = "/global/scratch/ekoch/combined/14B-088_channel_ms/"

start_chan = 11
nchan = 205

ms_split_by_channel(vis, nchan=nchan, start=start_chan,
                    output_dir=output_dir, datacolumn='DATA')
