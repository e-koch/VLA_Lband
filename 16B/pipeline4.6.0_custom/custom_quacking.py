
'''
Both 16B-242 and 16B-236 have essentially identical observing tracks.
There are continual bad integrations which appear to partially be the result of
not having all antennas on target (in addition to the quacking applied in the
pipeline).
'''

from tasks import flagdata, flagmanager

# All based on integration times of 2 sec


# The beginning of 3C48
flagdata(vis=ms_active, mode='quack', field="3C48", quackinterval=25,
         quackmode='beg', flagbackup=False)

# Beginning and end of each J0119+3210
flagdata(vis=ms_active, mode='quack', field="J0119+3210", quackinterval=6,
         quackmode='beg', flagbackup=False)
flagdata(vis=ms_active, mode='quack', field="J0119+3210", quackinterval=6,
         quackmode='endb', flagbackup=False)

# Now the pol cals, which always start with ~10 sec of crap.
flagdata(vis=ms_active, mode='quack', field="J0319+4130", quackinterval=10,
         quackmode='beg', flagbackup=False)

flagdata(vis=ms_active, mode='quack', field="3C138", quackinterval=10,
         quackmode='beg', flagbackup=False)

# And finally the sources
flagdata(vis=ms_active, mode='quack', intent="*TARGET*", quackinterval=5,
         quackmode='beg', flagbackup=False)
flagdata(vis=ms_active, mode='quack', intent="*TARGET*", quackinterval=5,
         quackmode='endb', flagbackup=False)

flagmanager(vis=ms_active, mode='save', versionname="custom_quacking",
            comment='Produced by custom_quacking.py')
