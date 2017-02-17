
'''
Flagging for 16B-242 taken on 10/14/16.

Scheduling block: 16B-242.sb32681213.eb32956712.57675.0472325463.continuum

No issues from operator's log.

'''

# BP scan has about a minute w/o all antennas on target.
default("flagdata")
flagdata(vis='16B-242.sb32681213.eb32956712.57675.0472325463.continuum.ms',
         mode='manual', scan='2', timerange="<01:12:40",
         flagbackup=False)

# Additional channels dominated by RFI
default("flagdata")
flagdata(vis='16B-242.sb32681213.eb32956712.57675.0472325463.continuum.ms',
         mode='manual', spw='2:81~82,',
         flagbackup=False)

# Save the custom flags
flagmanager(vis='16B-242.sb32681213.eb32956712.57675.0472325463.continuum.ms',
        mode='save',
        versionname="flags_from_16B-242.sb32681213.eb32949259.57673.06817918981.continuum.custom_flagging.py",
        comment='Flagging from script for systematic issues')
