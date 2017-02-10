
'''
Flagging for 16B-242 taken on 10/14/16.

Scheduling block: 16B-242.sb32681213.eb32956712.57675.0472325463.continuum

No issues from operator's log.

'''

from tasks import flagdata, flagmanager

vis_name = "16B-236.sb32623829.eb32998154.57703.975014375.continuum.ms"

# Log says bad things happened with ea04. The continuum reduction looks
# relatively fine, but I don't trust it (the lines are terrible).
# default("flagdata")
# flagdata(vis=vis_name,
#          mode='manual', antenna='ea04',
#          flagbackup=False)

# Individual scans, etc.
# Beginning of 3C48
default('flagdata')
flagdata(vis=vis_name, mode='manual', scan='2',
         timerange="<23:31:20", flagbackup=False)
# Scan 3 of the gain cal looks terrible
flagdata(vis=vis_name, mode='manual', scan='3',
         timerange="<23:34:50", flagbackup=False)

# Spike in scan 52
flagdata(vis=vis_name, mode='manual', scan='52',
         spw='2,5,7', timerange="2:32:08~2:33:10",
         flagbackup=False)

# Two spikes in SPW 2 that persist for a large portion of the track
# Remove first 40 channels to avoid it. Also remove these in 3C138
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='9~42,65',
         spw='2:0~41', flagbackup=False)
# A couple bad channels in scans 60,61 for spw 3
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='60~61',
         spw='3:119~', flagbackup=False)
# Longer bad spike in scan 52
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='52',
         spw='3', timerange="2:31:40~2:32:30",
         flagbackup=False)

# All sorts of time-variable crap in SPW 5
# One bad channel near known RFI for part of the track
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='5~44',
         spw='5:51', flagbackup=False)

# Scan 8 has extra RFI peak
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='8',
         spw='5:78~80', flagbackup=False)
# Scan 10 first channels crap
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='10',
         spw='5:0~21', flagbackup=False)
flagdata(vis=vis_name, mode='manual', scan='33',
         spw='5:71~74', flagbackup=False)
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='14',
         spw='5', timerange=">0:16:50",
         flagbackup=False)
flagdata(vis=vis_name, mode='manual', scan='54',
         spw='5', timerange="<2:39:00",
         flagbackup=False)
flagdata(vis=vis_name, mode='manual', scan='62',
         spw='5', timerange="03:10:15~03:10:30",
         flagbackup=False)
# There's a big spike that appears for only a few scans at a time.
default("flagdata")
flagdata(vis=vis_name, mode='manual',
         scan='10,15~16,36~38,62~64',
         spw='5:18~21',
         flagbackup=False)

# Shorter 52 spike in 6
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='52',
         spw='6', timerange="2:32:08~2:32:30",
         flagbackup=False)

# Large spike in scan 10 of spw 7
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='10',
         spw='7', timerange="<24:00:00",
         flagbackup=False)
flagdata(vis=vis_name, mode='manual', scan='10',
         spw='7', timerange="24:01:00~24:01:40",
         flagbackup=False)
# Flag ea21 in spw 7 J0319+4130 scan. Not sure what happened here...
default("flagdata")
flagdata(vis=vis_name, mode='manual', scan='64',
         spw='7', antenna='ea21',
         flagbackup=False)


# Save the custom flags
flagmanager(vis=vis_name,
            mode='save',
            versionname="manual_script_flagging",
            comment='Flagging from script for systematic issues')
