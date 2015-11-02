
'''
Extracts info from the CLEAN logs in separate channels.
'''

import sys
import glob
import os

from casa_tools import collect_clean_results


path_to_results = sys.argv[1]
output_filename = sys.argv[2] if sys.argv[2] != "None" else None
show_in_browser = True if sys.argv[3] == "True" else False

# In this case, each results is within it's own directory
# ie. path_to_results/channel_$i/
channel_direcs = glob.glob(os.path.join(path_to_results, "channel_*"))

# Now grab the logs
log_files = []
for channel in channel_direcs:
    log_file = glob.glob(os.path.join(channel, "casa*.log"))

    if not log_files:
        print("Cannot find log file in "+channel)
        continue

    # I'm going to assume that there should only be one logfile.
    # This should be guaranteed when using my automated setup.
    log_files.append(log_file[0])

# Now feed it into the results extractor
collect_clean_results(log_files, filename=output_filename,
                      show_in_browser=show_in_browser)
