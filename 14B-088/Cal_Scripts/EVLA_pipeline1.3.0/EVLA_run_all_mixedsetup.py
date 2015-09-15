
import sys
import os
import numpy as np
import copy
import glob
import shutil

'''
EVLA pipeline running for mixed setups.

Applies the first few steps of the pipeline up to basic flagging
(shadowing, zeros, ...). The continuum and line spws are then split
into their own directories. The full pipeline is then run on each. The lines
do not have rflag run on them.
'''

# Give the name of the ms. The folder containing the ms will be used for the
# splits
try:
    vis = sys.argv[1]
    path_to_pipeline = sys.argv[2]
    hanning_smooth = sys.argv[3]
except IndexError:
    vis = raw_input("MS File? : ")
    path_to_pipeline = raw_input("Path to pipeline? : ")
    hanning_smooth = raw_input("Hanning smooth? : ")

if vis[-1] == "/":
    vis = vis[:-1]

if not path_to_pipeline[-1] == "/":
    path_to_pipeline += "/"

# Chop off the .ms
SDM_name = vis[:-3]
SDM_name_orig = copy.copy(SDM_name)

# Set Hanning smoothing
if hanning_smooth == 'y':
    myHanning = 'y'
else:
    myHanning = 'n'

# Figure out which are the lines and which are the continuum SPWs.

tb.open(vis + '/SPECTRAL_WINDOW')
bandwidths = tb.getcol('TOTAL_BANDWIDTH')
tb.close()

tb.open(vis + '/FIELD')
fields = tb.getcol('NAME')
tb.close()

# Drop the pol cal.
fields = fields[np.where(fields != '0521+166=3C138')]

# Define a threshold between expected bandwidths
# Going with 10 MHz
thresh_bw = 1.0e7

spws = np.arange(0, len(bandwidths))

line_spws = [str(i) for i in spws[np.where(bandwidths < thresh_bw)]]
cont_spws = [str(i) for i in spws[np.where(bandwidths > thresh_bw)]]

print("Line SPWs: " + str(line_spws))
print("Coninuum SPWs: " + str(cont_spws))

print("Running initial pipeline.")

execfile(path_to_pipeline + "EVLA_pipeline_initial_mixed.py")

print("Splitting by SPW.")
os.mkdir('speclines')

split(vis=vis, outputvis="speclines/"+SDM_name+".speclines.ms",
      spw=",".join(line_spws), datacolumn='DATA', field=",".join(fields))

os.mkdir('continuum')

split(vis=vis, outputvis="continuum/"+SDM_name+".continuum.ms",
      spw=",".join(cont_spws), datacolumn='DATA', field=",".join(fields))

print("Running full pipeline on the spectral lines.")

os.chdir("speclines")

SDM_name = SDM_name_orig+".speclines"
myHanning = 'n'

execfile(path_to_pipeline + "EVLA_pipeline_lines.py")

# It saves a bunch of plots to the parent directory. So, move them to weblog here.
pngs = glob.glob("../*.png")
for png in pngs:
    shutil.move(png, "weblog/")

print("Running full pipeline on the spectral lines.")

os.chdir("../continuum")

SDM_name = SDM_name_orig+".continuum"
myHanning = 'n'

execfile(path_to_pipeline + "EVLA_pipeline_continuum.py")

pngs = glob.glob("../*.png")
for png in pngs:
    shutil.move(png, "weblog/")

print("All done!")
