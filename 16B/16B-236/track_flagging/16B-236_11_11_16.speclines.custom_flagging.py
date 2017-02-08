
'''
Flagging for 16B-242 taken on 10/14/16.

Scheduling block: 16B-242.sb32681213.eb32956712.57675.0472325463.continuum

No issues from operator's log.

'''

from tasks import flagdata, flagmanager

vis_name = "16B-236.sb32623829.eb32998154.57703.975014375.speclines.ms"

# Log says bad things happened with ea04. And is it ever right!
default("flagdata")
flagdata(vis=vis_name,
         mode='manual', antenna='ea04',
         flagbackup=False)

# Intermittent issues with ea15, ea17. Online flags seem to have
# removed this section already, but just to be sure.
# flagdata(vis=vis_name,
#          mode='manual', antenna='ea17',
#          timerange="23:24:00~03:23:40",
#          flagbackup=False)

flagdata(vis=vis_name,
         mode='manual', antenna='ea15',
         timerange="01:03:58~01:58:15",
         flagbackup=False)

# Individual scans, etc.
# Beginning of 3C48
default('flagdata')
flagdata(vis=vis_name, mode='manual', scan='2',
         timerange="<23:31:20", flagbackup=False)

# Bad baselines in scan 43? Find the antenna causing it
default('flagdata')
# flagdata(vis=vis_name, mode='manual', scan='43',
#          antenna='XXX', flagbackup=False)

default('flagdata')
# Big RFI peak over half of band
flagdata(vis=vis_name, mode='manual', scan='52',
         spw='0~1,7~9', timerange="2:32:08~2:32:30",
         flagbackup=False)
# Noise spike at beginning of track
flagdata(vis=vis_name, mode='manual', scan='10',
         spw='0,8', timerange='<23:50:00',
         flagbackup=False)
# Big noise spikes in 3 channels for scans 42,44,45
# flagdata(vis=vis_name, mode='manual', scan='42,44,45',
#          spw='0', channel="XXX",
#          flagbackup=False)
# flagdata(vis=vis_name, mode='manual', scan='42,44,45',
#          spw='1', channel="XXX",
#          flagbackup=False)
# Multiple time spikes in spw 3 (1612)
flagdata(vis=vis_name, mode='manual', scan='10',
         spw='3', timerange="<24:00:00",
         flagbackup=False)
flagdata(vis=vis_name, mode='manual', scan='12',
         spw='3', timerange=">00:08:30",
         flagbackup=False)
flagdata(vis=vis_name, mode='manual', scan='52',
         spw='3~6', timerange="2:32:08~2:33:10",
         flagbackup=False)
flagdata(vis=vis_name, mode='manual', scan='54',
         spw='3', timerange="<2:39:00",
         flagbackup=False)


# Save the custom flags
flagmanager(vis=vis_name,
            mode='save',
            versionname="manual_script_flagging",
            comment='Flagging from script for systematic issues')
