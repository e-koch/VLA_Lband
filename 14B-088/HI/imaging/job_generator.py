
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
#PBS -l pmem=2000m
#PBS -l feature=X5675
#PBS -l nodes=1:ppn=12
#PBS -l walltime=72:00:00
#PBS -l epilogue=/home/ekoch/code_repos/simscript/epilogue.sh

source /home/ekoch/.bashrc

cd X1

echo "Starting at: `date`"
casa --nologger --logfile X5 -c /home/ekoch/code_repos/VLA_Lband/14B-088/HI/imaging/HI_single_channel_clean.py X2 X3 X4
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


# Set the directory to look in.
ms_channel = "/home/ekoch/m33/14B-088/14B-088_channel_ms/"
model_channel_name = "/home/ekoch/m33/14B-088/model_channels/M33_14B-088_HI_" \
    "model_channel_{}.image"
mask_channel_name = "/home/ekoch/m33/14B-088/mask_channels/M33_14B-088_HI_" \
    "mask_channel_{}.image"
output_direc = "/home/ekoch/m33/14B-088/single_channels/"

# Use mask and model? Disable when continuing to clean.
use_mask_model = True

# Set the mode to use. Continuously checking for new splits, or a set number.
sub_mode = sys.argv[1]
if sub_mode == "continuous":
    pass
elif sub_mode == "range":
    try:
        start = int(sys.argv[2])
        stop = int(sys.argv[3])
    except IndexError:
        raise IndexError("Must provide a start and stop when using "
                         "'range' mode.")
else:
    raise TypeError("sub_mode must be 'continuous' or 'range'.")

while True:

    # Run channels in given range
    if sub_mode == "range":
        channel_ms = []
        for chan in xrange(start, stop + 1):
            channel_path = \
                os.path.join(ms_channel,
                             "14B-088_HI_LSRK.ms.contsub_channel_{}.ms".format(chan))
            channel_ms.append(channel_path)

    elif sub_mode == "continuous":
        channel_ms = glob.glob(os.path.join(ms_channel, "*channel*.ms"))

        channel_ms = drop_last(channel_ms)

    # If there aren't any more split ms in the path, break and exit
    if len(channel_ms) == 0:
        print("No more MSs found in the directory. Exiting.")
        break

    # Now loop through the existing channel ms
    for chan in channel_ms:
        chan_num = int(chan.split("_")[-1][:-3])
        # adjust for numbering offset
        mod_mask_num = chan_num - 670

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
        if not os.path.exists(model_name):
            shutil.move(model_channel_name.format(mod_mask_num), channel_direc)
        if not os.path.exists(mask_name):
            shutil.move(mask_channel_name.format(mod_mask_num), channel_direc)

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

    # Temporary stopper
    break

    # Wait an hour, then check again for new channel ms
    time.sleep(3600)
