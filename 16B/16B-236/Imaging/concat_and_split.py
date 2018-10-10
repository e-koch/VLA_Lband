
'''
Combine the tracks, then split out the science fields
'''

import os
from glob import glob

from tasks import virtualconcat, split

# Grab all of the MS tracks in the folder (should be 12)
myvis = glob("*.speclines.ms")

assert len(myvis) == 12

default('concat')
virtualconcat(vis=myvis, concatvis='16B-236_lines_all.ms',
              keepcopy=False)

default('split')
split(vis='16B-236_lines_all.ms', outputvis='16B-236_lines.ms',
      field='M33*',
      datacolumn='corrected',
      keepflags=False)

os.system("rm -r 17B-162_lines_all.ms")
