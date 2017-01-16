
'''
Template for adding in flagging for systematic/persistent issues in a track.

This example is for 16B-242.sb32681213.eb32949259.57673.06817918981.continuum
(16B-242_10_12_16).

Rename this script to {TRACK_NAME}.{continuum or speclines}.custom_flagging.py
and place in the pipeline directory. The script will be moved into the
restoration_products folder when the cleanup script is run.
'''

# ea21 had dropped DPM and returned at 5:02, but there still seems to be
# an issue with it in SPW 7. Flag ea21 after 5:02 in SPW 7.
default("flagdata")
flagdata(vis='16B-242.sb32681213.eb32949259.57673.06817918981.continuum.ms',
	 mode='manual', spw='7', antenna='ea21', timerange=">04:36:02",
	 flagbackup=False)

# Save the custom flags
flagmanager(vis='16B-242.sb32681213.eb32949259.57673.06817918981.continuum.ms',
	    mode='save',
        versionname="flags_from_16B-242.sb32681213.eb32949259.57673.06817918981.continuum.custom_flagging.py",
        comment='Flagging from script for systematic issues')
