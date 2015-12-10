
'''
Split the HI data from AT0206 into individual channels.
'''

from casa_tools import ms_split_by_channel

vis = "/global/scratch/ekoch/combined/M33_b_c_LSRK.ms"
output_dir = "/global/scratch/ekoch/combined/AT0206_channel_ms/"

start_chan = 10
nchan = 205

ms_split_by_channel(vis, nchan=nchan, start=start_chan,
                    output_dir=output_dir, datacolumn='DATA', spw='*')
