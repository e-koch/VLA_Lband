
'''
Split the HI Arecibo mask for 14B-088 into individual
channels using the undilated mask for certain channels.

These channels began to diverge with the new mask, but worked fine
with the old one.

Requires the results csv produced by channel_clean_results.csv
'''

from astropy.table import Table
import numpy as np

from casa_tools import image_split_by_channel


# Load the results csv file
t = Table.read("/home/eric/Dropbox/code_development/VLA_Lband/14B-088/HI/channel_clean_results.csv")

bad_chans = []

for i, row in enumerate(t):
    if row['Reached Threshold'] == "True":
        continue
    if np.isnan(row['Max Residual']):
        continue

    chan = int(row["Name"].split("/")[-2].split("_")[-1])
    bad_chans.append(chan-670)

print "Splitting now..."

mask = "M33_14B-088_HI_mask.image"
output_dir = "old_mask_channels/"

nchan = 1

for start_chan in bad_chans:
    print start_chan
    image_split_by_channel(mask, nchan=nchan, start=start_chan,
                           output_dir=output_dir)
