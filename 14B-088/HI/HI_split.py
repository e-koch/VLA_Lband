
'''
For splitting the HI data from 14B-088 data.
'''

import glob
import os

ms_path = raw_input("Path and name of MS?: ")
output_path = raw_input("Output name and path?: ")

spw = '1'
chans = '3000~5000'
uvsub_channels = '0:0~200;1860~2000'
field = 'M33*'

spw_slice = spw + ":" + chans

split(vis=ms_path, outputvis=output_path, spw=spw_slice,
      datacolumn='corrected', field=field,
      correlation='RR,LL', keepflags=False)

uvcontsub(vis=output_path, field=field,
          fitspw=uvsub_channels, solint='int',
          fitorder=0)
