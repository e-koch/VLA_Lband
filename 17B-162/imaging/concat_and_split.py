
'''
Combine the tracks, then split out the science fields
'''

import os
import sys
from glob import glob

from tasks import concat, split

# Grab all of the MS tracks in the folder (should be 17)
myvis = glob("*.speclines.ms")

assert len(myvis) == 17

default('concat')
concat(vis=myvis, concatvis='17B-162_lines_all.ms', timesort=False)

default('split')
split(vis='17B-162_lines_all.ms', outputvis='17B-162_lines.ms',
      field='M33*',
      datacolumn='corrected',
      keepflags=False)

# default('split')
# split(vis='17B-162_lines_all.ms',
#       outputvis='17B-162_lines_cals.ms',
#       field='J0319+4130,*3C138,*3C48',
#       datacolumn='corrected',
#       keepflags=False)

# os.system("rm -r 17B-162_lines_all.ms")
