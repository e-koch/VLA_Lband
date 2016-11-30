
'''
Generate a pbs script for job submission, submit the job, be happy
'''

import glob
import sys
import os
import shutil
import time


def return_template(output_direc, ms_name, model_name, mask_name):

    # Emailing
    #PBS -m bea
    #PBS -M koch.eric.w@gmail.com

    template = \
        '''
#!/bin/bash

#PBS -S /bin/bash
#PBS -l pmem=1000m
#PBS -l nodes=1:ppn=12
#PBS -l walltime=7:00:00
#PBS -l epilogue=/home/ekoch/code_repos/simscript/epilogue.sh

source /home/ekoch/.bashrc

cd X1

echo "Starting at: `date`"
casa-4.4 --nologger --logfile X5 -c /home/ekoch/code_repos/VLA_Lband/AT0206/Imaging/HI_single_channel_clean.py X2 X3 X4
echo "Exited with code $? at: `date`"
        '''

    template = template.strip()

    template = template.replace("X1", output_direc)
    template = template.replace("X2", ms_name)
    template = template.replace("X3", model_name)
    template = template.replace("X4", mask_name)

    # Create log file name
    logfile = ms_name.rstrip(".ms") + ".log"

    template = template.replace("X5", logfile)

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


# Set the path to the directory with the data
# This should contain the folders:
#       * channel_ms - folder with the channel MS splits
#       * model_channels - folder with split model channels
#       * mask_channels - same for the mask
data_directory = sys.argv[1]

# A common string input, based on the initial ms. This follows the output of
# the channel MS split and
file_prefix = sys.argv[2]

# Set the MS channel numbers to loop through (start, stop)
try:
    start = int(sys.argv[2])
    stop = int(sys.argv[3])
except IndexError:
    raise IndexError("Must provide a start and stop channel number.")

# Set the directory to look in.
ms_channel = os.path.join(data_directory, "channel_ms/")
model_channel_name = \
    os.path.join(data_directory,
                 "model_channels/M33_model_channel_{}.image")
mask_channel_name = \
    os.path.join(data_directory,
                 "mask_channels/M33_newmask_channel_{}.image")

output_direc = os.path.join(data_directory, "single_channels/")

# Use mask and model? Disable when continuing to clean.
use_mask = True
# Due to the current issues with the Arecibo data as a SD model, skip using it
use_model = False

# Run channels in given range
channel_ms = []
for chan in xrange(start, stop + 1):
    channel_path = \
        os.path.join(ms_channel,
                     "M33_b_c_LSRK_channel_{}.ms".format(chan))
    channel_ms.append(channel_path)

# If there aren't any more split ms in the path, break and exit
if len(channel_ms) == 0:
    print("No more MSs found in the directory. Exiting.")
    raise ValueError("No measurement sets found.")

# Now loop through the existing channel ms
for chan in channel_ms:
    chan_num = int(chan.split("_")[-1][:-3])
    # adjust for numbering offset
    mod_mask_num = chan_num - 10

    channel_direc = os.path.join(output_direc,
                                 "channel_{}".format(chan_num))

    # Check if that channel has been imaged already
    # if os.path.isdir(channel_direc):
    #     print("Already imaged "+str(chan_num)+". Skipping")
    #     continue

    if not os.path.exists(channel_direc):
        os.mkdir(channel_direc)

    # Names of the clean inputs in the channel_ms folder. Defined here to
    # check if they exist before moving.
    base_ms_name = os.path.basename(chan.rstrip("/"))
    chan_ms = os.path.join(channel_direc, base_ms_name)

    base_model_name = \
        os.path.basename(model_channel_name.format(mod_mask_num))
    model_name = os.path.join(channel_direc, base_model_name)
    base_mask_name = \
        os.path.basename(mask_channel_name.format(mod_mask_num))
    mask_name = os.path.join(channel_direc, base_mask_name)

    # Now move the mask, model, and channel ms into the folder, if they
    # aren't there already
    if not os.path.exists(chan_ms):
        shutil.move(chan, channel_direc)
    if not os.path.exists(model_name) and use_model:
        shutil.move(model_channel_name.format(mod_mask_num), channel_direc)
    if not os.path.exists(mask_name) and use_mask:
        shutil.move(mask_channel_name.format(mod_mask_num), channel_direc)

    if not use_mask:
        mask_name = "None"
    if not use_model:
        model_name = "None"

    chan_template = return_template(channel_direc, chan_ms,
                                    model_name, mask_name)

    # Write to file
    sub_file = os.path.join(channel_direc,
                            "channel_{}.sub".format(chan_num))
    if not os.path.exists(sub_file):
        with open(sub_file, 'w') as f:
            f.write(chan_template)

    # Now submit!
    old_direc = os.getcwd()
    os.chdir(channel_direc)  # Switch to directory so log files are there
    os.system("qsub " + sub_file)
    os.chdir(old_direc)
