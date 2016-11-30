
'''
Generate a pbs script for job submission, submit the job, be happy
'''

import glob
import sys
import os
import shutil
import time


def return_template(output_direc, ms_name_1, ms_name_2, model_name,
                    mask_name, outname):

    # Emailing
    # PBS -m bea
    # PBS -M koch.eric.w@gmail.com

    template = \
        '''
#!/bin/bash

#PBS -S /bin/bash
#PBS -l pmem=1000m
#PBS -l nodes=1:ppn=12
#PBS -l walltime=7:00:00

source /home/ekoch/.bashrc

cd X1

echo "Starting at: `date`"
casa-4.4 --nologger --logfile X7 -c /home/ekoch/VLA_Lband/14B-088/HI/archival_combination/HI_single_channel_clean.py X2 X3 X4 X5 X6
echo "Exited with code $? at: `date`"
        '''

    template = template.strip()

    template = template.replace("X1", output_direc)
    template = template.replace("X2", ms_name_1)
    template = template.replace("X3", ms_name_2)
    template = template.replace("X4", model_name)
    template = template.replace("X5", mask_name)
    template = template.replace("X6", outname)

    # Create log file name
    logfile = outname + ".log"

    template = template.replace("X7", logfile)

    return template


def drop_last(ms_list):
    '''
    CASA is generally writing to the final MS in the folder, so skip it.
    '''

    max_num = 0
    for ms in ms_list:
        if int(ms.split("_")[-1][:-3]) > max_num:
            max_num_ms = ms
            max_num = int(ms.split("_")[-1][:-3])

    ms_list.remove(max_num_ms)

    return ms_list


# Set the directory to look in.
data_direc = sys.argv[1]

ms_channel_1 = os.path.join(data_direc, "14B-088_channel_ms/")
ms_channel_2 = os.path.join(data_direc, "AT0206_channel_ms/")
model_channels = os.path.join(
    data_direc, "model_channels/M33_model_channel_")
mask_channels = os.path.join(
    data_direc, "mask_channels/M33_newmask_channel_")
output_direc = os.path.join(data_direc, "single_channels/")

# Set the mode to use. Continuously checking for new splits, or a set number.
try:
    start = int(sys.argv[2])
    stop = int(sys.argv[3])
except IndexError:
    raise IndexError("Must provide a start and stop channel number.")

# Use mask and model? Disable when continuing to clean.
use_mask = True
# Due to the current issues with the Arecibo data as a SD model, skip using it
use_model = False

# Run channels in given range
channel_ms_1 = []
channel_ms_2 = []
for chan in xrange(start, stop + 1):
    channel_path = \
        os.path.join(ms_channel_1,
                     "14B-088_HI_LSRK_AT0206_regrid.ms.contsub_channel_" +
                     str(chan) + ".ms")
    channel_ms_1.append(channel_path)

    channel_path = \
        os.path.join(ms_channel_2,
                     "M33_b_c_LSRK_channel_" +
                     str(chan) + ".ms")
    channel_ms_2.append(channel_path)

# If there aren't any more split ms in the path, break and exit
if len(channel_ms_1) == 0 or len(channel_ms_2) == 0:
    raise ValueError("No more MSs found in the directory. Exiting.")


# Now loop through the existing channel ms
for chan_1, chan_2 in zip(channel_ms_1, channel_ms_2):
    chan_num_1 = int(chan_1.split("_")[-1][:-3])
    chan_num_2 = int(chan_2.split("_")[-1][:-3])

    assert chan_num_1 == chan_num_2
    chan_num = chan_num_1

    # adjust for numbering offset
    mod_mask_num = chan_num - 10

    channel_direc = os.path.join(output_direc, "channel_" + str(chan_num))

    if not os.path.exists(channel_direc):
        os.mkdir(channel_direc)
    if not os.path.exists(chan_1):
        shutil.move(chan_1, channel_direc)
    if not os.path.exists(chan_2):
        shutil.move(chan_2, channel_direc)

    model_name = model_channels + str(mod_mask_num) + ".image"
    if not os.path.exists(channel_direc) and use_model:
        shutil.move(model_name, channel_direc)

    mask_name = mask_channels + str(mod_mask_num) + ".image"
    if not os.path.exists(channel_direc) and use_mask:
        shutil.move(mask_name, channel_direc)

    chan_ms_1 = os.path.join(channel_direc, chan_1.split("/")[-1])
    chan_ms_2 = os.path.join(channel_direc, chan_2.split("/")[-1])

    if use_model:
        model_name = os.path.join(channel_direc,
                                  "M33_model_channel_"
                                  + str(mod_mask_num) + ".image")
    else:
        model_name = "None"

    if use_mask:
        mask_name = os.path.join(channel_direc,
                                 "M33_newmask_channel_"
                                 + str(mod_mask_num) + ".image")
    else:
        mask_name = "None"

    outname = \
        os.path.join(channel_direc,
                     "M33_HI_14B-088_AT0206_channel_" + str(chan_num))

    chan_template = return_template(channel_direc, chan_ms_1, chan_ms_2,
                                    model_name, mask_name, outname)

    # Write to file
    sub_file = os.path.join(channel_direc, "channel_" + str(chan_num) + ".sub")
    if not os.path.exists(sub_file):
        with open(sub_file, 'w') as f:
            f.write(chan_template)

    # Now submit!
    old_direc = os.getcwd()
    os.chdir(channel_direc)  # Switch to directory so log files are there
    os.system("qsub " + sub_file)
    os.chdir(old_direc)
