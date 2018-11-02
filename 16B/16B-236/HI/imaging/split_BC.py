
'''
Split the 14B and 17B fields that overlap with 16B-236 pointing.
'''

import os

from tasks import split

fields = ['M33_7_center', 'M33_6', 'M33_14', 'M33_12']

fourteenB_ms = "14B-088_HI_LSRK.ms.contsub"
fourteenB_ms_out = "14B-088_HI_LSRK.16B236_fields.ms.contsub"

seventeenB_ms = "17B-162_HI_spw_0_LSRK.ms.contsub"
seventeenB_ms_out = "17B-162_HI_spw_0_LSRK.16B236_fields.ms.contsub"

data_path = os.path.expanduser("~/space/ekoch/VLA_tracks/17B-162/products/HI/")

split(vis=os.path.join(data_path, fourteenB_ms),
      outputvis=fourteenB_ms_out,
      datacolumn='data',
      field=",".join(fields))

split(vis=os.path.join(data_path, seventeenB_ms),
      outputvis=seventeenB_ms_out,
      datacolumn='data',
      field=",".join(fields))
