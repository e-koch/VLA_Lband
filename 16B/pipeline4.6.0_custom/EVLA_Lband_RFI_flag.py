
'''
Several portions of the L-band have constant RFI that is unrecoverable.

This is meant to be run in the pipeline environment to remove these regions
before wasting time trying to find solutions.
'''

from glob import glob
import os

from tasks import flagdata, flagmanager


logprint("Starting EVLA_Lband_RFI_flag.py", logfileout='logs/known_RFI_flag.log')

try:
    ms_active
except NameError:
    raise NameError("Run EVLA_pipe_restore.py")

# Check if the flag version already exists. If it does, skip the flagging and
# continue
needs_flagging = True
flag_folder = "{}.flagversions".format(ms_active)
if os.path.exists(flag_folder):
    flag_versions = \
        glob(os.path.join(flag_folder, "flags.known_RFI*"))
    if len(flag_versions) != 0:
        needs_flagging = False

if needs_flagging:
    # Flag SPW 0, 1 and 4. The pipeline always fails for 0 and usually for 4.
    # SPW 1 is always returned by the pipeline, but by-eye investigations and
    # the VLA RFI plots show there are almost no RFI free channels.
    default("flagdata")
    flagdata(vis=ms_active, spw="0,1,4", flagbackup=False)

    # Now remove specific portions of the other SPWs
    # Channel ranges are for 16B, and are the regions that were found to
    # persist over multiple SBs. Other strong RFI exists (SPW 6, for instance)
    # but is intermittent.
    # This leaves 438 good channels, or 68% of the recoverable 5 SPWs.
    # This is 43% of the total BW. We assumed 60% recovered in the proposal.
    # Using the sensitivity calculator, the difference in sensitivity is
    # 1.88 uJy/bm vs. 1.63 uJy/bm over 48 hours. So not much of a difference.
    spw2 = "2:0~20;42~54;83~96"
    spw3 = "3:0~10;30;31;48;49;69;70"
    spw5 = "5:52~67;112;113"
    # spw6 = ""  # There are two narrow, strong, but intermittent RFI sources
    spw7 = "7:44~"

    flag_str = ",".join([spw1, spw2, spw3, spw5, spw7])
    flagdata(vis=ms_active, spw=flag_str, flagbackup=False)

    flagmanager(vis=ms_active, mode='save', versionname="known_RFI")
else:
    logprint("known_RFI flag version already exists. Skipping flagging.")

logprint("Finished EVLA_Lband_RFI_flag.py", logfileout='logs/known_RFI_flag.log')
