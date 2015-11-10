
'''
Convert each MS to an mms
'''

from tasks import partition

vis = "/home/ekoch/scratch/combined/14B-088_HI_LSRK_AT0206_regrid.ms.contsub"
outputvis = \
    "/home/ekoch/scratch/combined/14B-088_HI_LSRK_AT0206_regrid.ms.contsub.mms"

partition(vis=vis, outputvis=outputvis)


vis = "/home/ekoch/scratch/combined/M33_b_c_LSRK.ms"
outputvis = \
    "/home/ekoch/scratch/combined/M33_b_c_LSRK.mms"

partition(vis=vis, outputvis=outputvis)
