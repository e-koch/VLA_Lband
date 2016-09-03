
'''
Remove calibration and statwt from an MS. Optionally remove all flagging.
'''

import os
import sys

ms_name = raw_input("MS Name: ")
remove_flagging = \
    True if raw_input("Remove all flagging (y/n): ") == 'y' else False

print("Removing calibration. This may take awhile.")
default("clearcal")
clearcal(vis=ms_name, addmodel=False)

print("Clearing statwts")
default("clearstat")
clearstat(vis=ms_name)

if remove_flagging:
    print("Removing all flags.")
    default("flagdata")
    flagdata(vis=ms_name, mode='unflag')

print("Finished resetting MS. Note that that data are still Hanning smoothed,"
      " if this pipeline option was chosen.")
