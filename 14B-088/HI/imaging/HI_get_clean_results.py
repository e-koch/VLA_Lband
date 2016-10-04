
'''
Extracts info from the CLEAN logs in separate channels.
'''

import sys
import glob
import os
import numpy as np

# from casa_tools import collect_clean_results
execfile("/home/eric/Dropbox/code_development/VLA_Lband/CASA_tools/casa_tools/extract_from_log.py")


path_to_results = sys.argv[1]
output_filename = sys.argv[2] if sys.argv[2] != "None" else None
log_prefix = sys.argv[3]  # The common log file starting string
show_in_browser = True if sys.argv[4] == "True" else False

# In this case, each results is within it's own directory
# ie. path_to_results/channel_$i/
channel_direcs = glob.glob(os.path.join(path_to_results, "channel_*"))

# Sort by the channel number
min_chan = 10000
max_chan = 0

for channel in channel_direcs:
    # Last 3 characters
    chan_num = int(channel.split("_")[-1])
    if chan_num < min_chan:
        min_chan = chan_num

    if chan_num > max_chan:
        max_chan = chan_num

ordered_chans = []
for chan in xrange(min_chan, max_chan+1):
    i = 0
    while True:
        chan_num = int(channel_direcs[i].split("_")[-1])
        if chan_num == chan:
            ordered_chans.append(channel_direcs[i])
            break
        i += 1
    else:
        print("No log found for channel "+str(chan))

channel_direcs = ordered_chans

# Now grab the logs
log_files = []
for channel in channel_direcs:
    log_file = glob.glob(os.path.join(channel, "{}*.log".format(log_prefix)))
    if not log_file:
        print("Cannot find log file in "+channel)
        continue

    # I'm going to assume that there should only be one logfile.
    # This should be guaranteed when using my automated setup.
    log_files.append(log_file[0])

# Now feed it into the results extractor
collect_clean_results(log_files, filename=output_filename,
                      show_in_browser=show_in_browser)
