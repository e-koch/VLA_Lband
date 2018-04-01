
import sys
import os
from tasks import split, importasdm

'''
Identify the continuum and line SPWs and split into separate MSs and
directories
'''

mySDM = sys.argv[-2]
parallel_run = True if sys.argv[-2] == "T" else False

ms_active = mySDM + ".ms"

importasdm(asdm=mySDM, vis=ms_active, ocorr_mode='co',
           applyflags=True, savecmds=True, tbuff=1.5,
           outfile='{}.flagonline.txt'.format(mySDM),
           createmms=parallel_run)

line_run = False
cont_run = False

parentdir = os.getcwd().split("/")[-1]

lines_folder = parentdir + '_speclines'
if not os.path.exists(lines_folder):
    os.mkdir(lines_folder)

split(vis=ms_active,
      outputvis=lines_folder + "/" + mySDM + ".speclines.ms",
      spw="8~17", datacolumn='DATA',
      field="0137+331=3C48,M33_2,M33_14,M33_6,M33_7_center,M33_12,M33_11,M33_8")

line_run = True

line_ms = mySDM + ".speclines"


cont_folder = parentdir + '_continuum'
if not os.path.exists(cont_folder):
    os.mkdir(cont_folder)

split(vis=ms_active,
      outputvis=cont_folder + "/" + mySDM + ".continuum.ms",
      spw="0~7", datacolumn='DATA',
      field="")

cont_ms = mySDM + ".continuum"
