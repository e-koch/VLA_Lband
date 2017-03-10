
'''
Due to a malfunction, this track was split in two. After the initial flagging,
they need to be combined before running the rest of the reduction.

The pre-amble is from EVLA_pipeline.py. See the attached license and version
changes there.
'''

import sys
import os
import time
import shutil

# Import path to pipeline
from paths import sixteenB_pipe_path
# Define location of pipeline
pipepath = sixteenB_pipe_path

version = "1.3.10"
svnrevision = '11nno'
date = "2017Feb10"

print "Pipeline version " + version + " for use with CASA 4.7.1"

[major, minor, revision] = casadef.casa_version.split('.')
casa_version = 100 * int(major) + 10 * int(minor) + int(revision)
if casa_version < 471:
    sys.exit("Your CASA version is " + casadef.casa_version +
             ", please re-start using CASA 4.7.1")
if casa_version > 471:
    sys.exit("Your CASA version is " + casadef.casa_version +
             ", please re-start using CASA 4.7.1")


# This is the default time-stamped casa log file, in case we
#    need to return to it at any point in the script
log_dir = 'logs'
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

maincasalog = casalogger.func_globals['thelogfile']


def logprint(msg, logfileout=maincasalog):
    print(msg)
    casalog.setlogfile(logfileout)
    casalog.post(msg)
    casalog.setlogfile(maincasalog)
    casalog.post(msg)
    return

# Create timing profile list and file if they don't already exist
if 'time_list' not in globals():
    time_list = []

timing_file = 'logs/timing.log'

if not os.path.exists(timing_file):
    timelog = open(timing_file, 'w')
else:
    timelog = open(timing_file, 'a')


def runtiming(pipestate, status):
    '''Determine profile for a given state/stage of the pipeline
    '''
    time_list.append(
        {'pipestate': pipestate, 'time': time.time(), 'status': status})
#
    if (status == "end"):
        timelog = open(timing_file, 'a')
        timelog.write(
            pipestate + ': ' + str(time_list[-1]['time'] -
                                   time_list[-2]['time']) + ' sec \n')
        timelog.flush()
        timelog.close()
        # with open(maincasalog, 'a') as casalogfile:
        #    tempfile = open('logs/'+pipestate+'.log','r')
        #    casalogfile.write(tempfile.read())
        #    tempfile.close()
        # casalogfile.close()
#
    return time_list

######################################################################

# The following script includes all the definitions and functions and
# prior inputs needed by a run of the pipeline.

SDM_name = "16B-242.sb32614458.eb32984320.57697.291263148145"

time_list = runtiming('startup', 'start')
execfile(pipepath + 'EVLA_pipe_startup.py')
time_list = runtiming('startup', 'end')
pipeline_save()

# IMPORT THE DATA TO CASA

execfile(pipepath + 'EVLA_pipe_import.py')

# GET SOME INFORMATION FROM THE MS THAT WILL BE NEEDED LATER, LIST
# THE DATA, AND MAKE SOME PLOTS

execfile(pipepath + 'EVLA_pipe_msinfo.py')

# DETERMINISTIC FLAGGING:
# TIME-BASED: online flags, shadowed data, zeroes, pointing scans, quacking
# CHANNEL-BASED: end 5% of channels of each spw, 10 end channels at
# edges of basebands

execfile(pipepath + 'EVLA_pipe_flagall.py')

# Copy the logs folder to a unique name
shutil("logs", "logs_1")

shutil.move("*.png", "logs_1")
shutil.move("*.listobs", "logs_1")

# Now the second part

SDM_name = "16B-242_sb32614458_5.57697.30203707176"

execfile(pipepath + 'EVLA_pipe_import.py')

execfile(pipepath + 'EVLA_pipe_msinfo.py')

execfile(pipepath + 'EVLA_pipe_flagall.py')

# Copy the logs folder to a unique name
shutil("logs", "logs_2")

shutil.move("*.png", "logs_2")
shutil.move("*.listobs", "logs_2")


# Combine the output MS's
concat(vis=["16B-242.sb32614458.eb32984320.57697.291263148145.ms",
            "16B-242_sb32614458_5.57697.30203707176.ms"],
       concatvis="16B-242.sb32614458.eb32984320.57697.291263148145.ms",
       timesort=True)
default("concat")

# Now run the pipeline as normal
SDM_name = "16B-242.sb32614458.eb32984320.57697.291263148145.ms"

try:

    ######################################################################

    # IMPORT THE DATA TO CASA

    execfile(pipepath + 'EVLA_pipe_import.py')

    # GET SOME INFORMATION FROM THE MS THAT WILL BE NEEDED LATER, LIST
    # THE DATA, AND MAKE SOME PLOTS

    execfile(pipepath + 'EVLA_pipe_msinfo.py')

    # DETERMINISTIC FLAGGING:
    # TIME-BASED: online flags, shadowed data, zeroes, pointing scans, quacking
    # CHANNEL-BASED: end 5% of channels of each spw, 10 end channels at
    # edges of basebands

    execfile(pipepath + 'EVLA_pipe_flagall.py')

    # Split into continuum and lines

    execfile(pipepath + 'EVLA_pipe_mixed_setup_split.py')

    current_path = os.getcwd()

    if line_run:
        line_path = line_settings["folder_name"]
        if not os.path.exists(line_path):
            raise Warning("speclines path was not created.")

        # Copy the online flags plot into each split directory
        os.system("cp ../onlineFlags.png .")

        os.chdir(line_path)

        # Change the SDM name and other needed settings
        SDM_name = line_settings['SDM_Name']
        myHanning = line_settings['myHanning']
        scratch = line_settings['scratch']

        execfile(pipepath + 'EVLA_pipeline_lines.py')

    os.chdir(current_path)

    if cont_run:
        cont_path = cont_settings["folder_name"]
        if not os.path.exists(cont_path):
            raise Warning("continuum path was not created.")

        # Copy the online flags plot into each split directory
        os.system("cp ../onlineFlags.png .")

        os.chdir(cont_path)

        # Change the SDM name and other needed settings
        SDM_name = cont_settings['SDM_Name']
        myHanning = cont_settings['myHanning']
        scratch = cont_settings['scratch']

        execfile(pipepath + 'EVLA_pipeline_continuum.py')

# Quit if there have been any exceptions caught:

except KeyboardInterrupt, keyboardException:
    logprint("Keyboard Interrupt: " + str(keyboardException))
except Exception, generalException:
    logprint("Exiting script: " + str(generalException))
