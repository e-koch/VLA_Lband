
'''
On mixed setups, flagall is already run. This defines the variables set
during that script so it doesn't need to be run multiple times.

BUT, if we're re-running the pipeline, the script will check for manual
flagging scripts saved in the appropriate project folder for the track.

'''

import os
import glob
import numpy as np

from paths import a_path

logprint("Starting EVLA_pipe_fake_flagall.py",
         logfileout='logs/flagall.log')
time_list = runtiming('flagall', 'start')
QA2_flagall = 'Pass'

logprint("These value are based on the current flags on the MS. All"
         " deterministic was performed prior to splitting into separate line"
         " and continuum MSs.")


# Path to the manual flagging scripts
# Do any of them match?
folder_name = os.getcwd().split("/")[-1]
proj_code = "16B-242" if "16B-242" in folder_name else "16B-236"
flagging_path = os.path.join(a_path, proj_code, "track_flagging")
track_flags = glob.glob(flagging_path + "/*.py")
hits = np.where([os.getcwd().split("/")[-1] in track for
                 track in track_flags])[0]
if len(hits) == 1:
    logprint("Found a manual flagging script to run: "
             "{}".format(hits[0].split(",")[-1]))
    execfile(track_flags[hits[0]])
elif len(hits) > 1:
    from warnings import warn
    import sys
    logprint("Multiple script hits? Something is wrong!")
    sys.exit()
else:
    logprint("No manual flagging script found.")


# report initial statistics
default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
myinitialflags = flagdata()
#clearstat()
logprint ("Initial flags summary", logfileout='logs/flagall.log')

start_total = myinitialflags['total']
start_flagged = myinitialflags['flagged']
logprint("Initial flagged fraction = {}".format(start_flagged / start_total),
         logfileout='logs/flagall.log')


init_on_source_vis = start_total

afterzero_total = start_total
afterzero_flagged = start_flagged

zero_flagged = 0.0

aftershadow_total = start_total
aftershadow_flagged = start_flagged

shadow_flagged = 0.0

flagdata_list = []
cmdreason_list = []

frac_flagged_on_source1 = 1.0 - ((start_total - start_flagged)  / init_on_source_vis)

if (frac_flagged_on_source1 >= 0.3):
    logprint("Initial flagging appears extensive. Check pipeline results to "
             "ensure usable solution are attained.")
    QA2_flagall='Fail'

logprint("Finished EVLA_pipe_fake_flagall.py", logfileout='logs/flagall.log')
logprint("QA2 score: "+QA2_flagall, logfileout='logs/flagall.log')
time_list = runtiming('flagall', 'end')

pipeline_save()
