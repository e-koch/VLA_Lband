
'''
Split the HI Arecibo mask for 14B-088 into individual
channels using the undilated mask for certain channels.

These channels began to diverge with the new mask, but worked fine
with the old one.

Requires the results csv produced by channel_clean_results.csv
'''

import numpy as np
from astropy.table import Table
import os

from casa_tools import image_split_by_channel
from tasks import rmtables


# Load the results csv file
t = Table.read("/home/eric/Dropbox/code_development/VLA_Lband/14B-088/HI/imaging/channel_clean_results_1016_4.4.csv")

# Do you want to delete the image products on Jasper, and copy the changed
# mask?
replace_mask_and_clean = True

bad_chans = []
bad_chans_ms = []

for i, row in enumerate(t):
    if row['Reached Threshold'] == "True":
        continue
    if np.isnan(row['Max Residual']):
        continue

    chan = int(row["Name"].split("/")[-2].split("_")[-1])
    bad_chans.append(chan - 670)
    bad_chans_ms.append(chan)

print "Splitting now..."

mask_path = os.path.expanduser("~/Data_3/M33/Arecibo/14B-088_items/")
mask = os.path.join(mask_path, "M33_14B-088_HI_mask.image")
output_dir = os.path.join(mask_path, "old_mask_channels/")

nchan = 1

for start_chan in bad_chans:
    print start_chan
    image_split_by_channel(mask, nchan=nchan, start=start_chan,
                           output_dir=output_dir)

np.save("bad_chans.npy", np.array(bad_chans_ms))

if replace_mask_and_clean:

    import shutil

    if not os.path.exists(os.path.expanduser("~/Jasper/m33")):
        raise Exception("Jasper doesn't seem to be mounted.")

    data_path = os.path.expanduser("~/Jasper/m33/14B-088/single_channels")

    for chan, chan_mask in zip(bad_chans_ms, bad_chans):
        print("Channel {}".format(chan))
        channel_path = os.path.join(data_path, "channel_{}".format(chan))

        # Remove old clean outputs
        rmtables(os.path.join(channel_path,
                              "14B-088_HI_LSRK.ms.contsub_channel_{}.clean.*".format(chan)))

        # log files, job outputs
        os.system("rm {0}/*.log {0}/*.last, {0}/*.sub.*".format(channel_path))

        # Now delete the old mask and copy over the new one.
        mask_name = os.path.join(channel_path,
                                 "M33_14B-088_HI_mask_channel_{}.image".format(chan_mask))
        rmtables(mask_name)

        new_mask_name = os.path.join(output_dir, os.path.basename(mask_name))
        shutil.copytree(new_mask_name, mask_name)
