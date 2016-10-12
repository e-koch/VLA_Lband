
import numpy as np
import os
from copy import copy
from taskinit import tb
from tasks import split


'''
Identify the continuum and line SPWs and split into separate MSs and
directories
'''

logprint("Starting EVLA_pipe_mixed_setup_split.py",
         logfileout='logs/mixed_setup_split.log')

# Figure out which are the lines and which are the continuum SPWs.
tb.open(ms_active + '/SPECTRAL_WINDOW')
bandwidths = tb.getcol('TOTAL_BANDWIDTH')
tb.close()

tb.open(ms_active + '/FIELD')
fields = tb.getcol('NAME')
tb.close()

# Define a threshold between expected bandwidths
# Going with 10 MHz
thresh_bw = 1.0e7

spws = np.arange(0, len(bandwidths))

line_spws = [str(i) for i in spws[np.where(bandwidths < thresh_bw)]]
cont_spws = [str(i) for i in spws[np.where(bandwidths > thresh_bw)]]

logprint("Line SPWS: {}".format(line_spws),
         logfileout='logs/mixed_setup_split.log')

logprint("Continuum SPWS: {}".format(cont_spws),
         logfileout='logs/mixed_setup_split.log')

logprint("Splitting by SPW.",
         logfileout='logs/mixed_setup_split.log')

line_run = False
cont_run = False

parentdir = os.getcwd().split("/")[-1]

if len(line_spws) == 0:
    logprint("No line SPWs found.",
             logfileout='logs/mixed_setup_split.log')
else:
    lines_folder = parentdir + '_speclines'
    if not os.path.exists(lines_folder):
        os.mkdir(lines_folder)

    line_fields = copy(fields)
    line_fields = line_fields[line_fields != "3C138"]
    line_fields = line_fields[line_fields != "J0319+4130"]

    split(vis=ms_active,
          outputvis=lines_folder + "/" + SDM_name + ".speclines.ms",
          spw=",".join(line_spws), datacolumn='DATA', field=",".join(fields))

    line_run = True

    line_settings = {"SDM_Name": SDM_name + ".speclines",
                     "scratch": False,
                     "myHanning": 'n',
                     "folder_name": lines_folder}

if len(cont_spws) == 0:
    logprint("No continuum SPWs found.",
             logfileout='logs/mixed_setup_split.log')
else:
    cont_folder = parentdir + '_continuum'
    if not os.path.exists(cont_folder):
        os.mkdir(cont_folder)

    split(vis=ms_active,
          outputvis=cont_folder + "/" + SDM_name + ".continuum.ms",
          spw=",".join(cont_spws), datacolumn='DATA', field=",".join(fields))

    cont_run = True

    cont_settings = {"SDM_Name": SDM_name + ".continuum",
                     "scratch": False,
                     "myHanning": myHanning,
                     "folder_name": cont_folder}

logprint("Finished EVLA_pipe_mixed_setup_split.py",
         logfileout='logs/mixed_setup_split.log')
