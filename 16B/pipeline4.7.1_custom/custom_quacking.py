
'''
Both 16B-242 and 16B-236 have essentially identical observing tracks.
There are continual bad integrations which appear to partially be the result of
not having all antennas on target (in addition to the quacking applied in the
pipeline).
'''

import numpy as np
import os

from tasks import flagdata, flagmanager

from paths import root

# All based on integration times of 2 sec
logprint("Beginning custom quacking")

default("flagdata")

# The beginning of 3C48
# This varies by when the track was started and all antennas got on source.
# Load in csv files with the start times.
if ms_active[:7] == "16B-242":
    flag_script_path = os.path.join(root, "16B/16B-242/track_flagging/")
    start_file = os.path.join(flag_script_path, "initial_on_source_times.csv")
    source_name = "NGC604"
elif ms_active[:7] == "16B-236":
    flag_script_path = os.path.join(root, "16B/16B-236/track_flagging/")
    start_file = os.path.join(flag_script_path, "initial_on_source_times.csv")
    source_name = "M33_Sarm"
else:
    raise ValueError("ms_active must use the default name starting with the"
                     " project ID (16B-242 or 16B-236).")

# Open the file and search for the track name
output = np.loadtxt(start_file, delimiter=',', dtype=str)

# Search for the track name
pos = np.where(output[:, 1] == ms_active[:-3])[0]

if len(pos) == 0:
    logprint("Could not find track name! Skipping 3C48 quack ...")
else:
    start_time = output[pos][0][-1]

    flagdata(vis=ms_active, mode='manual', field="3C48",
             flagbackup=False, scan='2', timerange="<{}".format(start_time))

# Beginning and end of each J0119+3210
# Adjust for slew time at beginning
flagdata(vis=ms_active, mode='quack', field="J0119+3210", quackinterval=13,
         quackmode='beg', flagbackup=False)
flagdata(vis=ms_active, mode='quack', field="J0119+3210", quackinterval=6,
         quackmode='endb', flagbackup=False)

# Now the pol cals, which always start with ~10 sec of crap. But we need to
# account for the slew time, which is the root of this issue: the online flags
# don't seem to go on long enough.
# For J0319+4130, the slew time is around 1:10 min
flagdata(vis=ms_active, mode='quack', field="J0319+4130", quackinterval=70,
         quackmode='beg', flagbackup=False)

# It's usually a longer slew to 3C138. Go with 1:45 min
flagdata(vis=ms_active, mode='quack', field="3C138", quackinterval=105,
         quackmode='beg', flagbackup=False)

# And finally the sources
# Adjust for slew time at beginning
flagdata(vis=ms_active, mode='quack', field=source_name, quackinterval=15,
         quackmode='beg', flagbackup=False)
flagdata(vis=ms_active, mode='quack', field=source_name, quackinterval=5,
         quackmode='endb', flagbackup=False)

flagmanager(vis=ms_active, mode='save', versionname="custom_quacking",
            comment='Produced by custom_quacking.py')

logprint("End custom quacking")
