
'''
Combine the tracks, then split out the science fields
'''

import os
from glob import glob

from tasks import concat, split

# Grab all of the MS tracks in the folder (should be 12)
myvis = glob("*.speclines.ms")

assert len(myvis) == 12

default('concat')
concat(vis=myvis, concatvis='16B-236_lines_all.ms')

default('split')
split(vis='16B-236_lines_all.ms', outputvis='16B-236_lines.ms',
      field='M33*',
      datacolumn='corrected',
      keepflags=False)

os.system("rm -r 16B-236_lines_all.ms")
