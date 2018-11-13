
'''
Combine the tracks, then split out the science fields
'''

import os
from glob import glob

from tasks import virtualconcat, split

# Grab all of the MS tracks in the folder (should be 12)
myvis = glob("16B-242.*.ms")

assert len(myvis) == 12

default('virtualconcat')
virtualconcat(vis=myvis, concatvis='16B-242_lines_all.ms',
              keepcopy=False)

default('split')
split(vis='16B-242_lines_all.ms', outputvis='16B-242_lines.ms',
      field='M33*',
      datacolumn='corrected',
      keepflags=False)

os.system("rm -r 16B-242_lines_all.ms")
