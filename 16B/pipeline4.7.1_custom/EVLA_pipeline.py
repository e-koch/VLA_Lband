######################################################################
#
# Copyright (C) 2013
# Associated Universities, Inc. Washington DC, USA,
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning VLA Pipelines should be addressed as follows:
#    Please register and submit helpdesk tickets via: https://help.nrao.edu
#    Postal address:
#              National Radio Astronomy Observatory
#              VLA Pipeline Support Office
#              PO Box O
#              Socorro, NM,  USA
#
######################################################################
# EVLA pipeline
# For continuum modes (contiguous spws within a baseband)
# May work for other modes as well
#
# 06/13/12 C. Chandler
# 07/20/12 B. Kent
# 02/05/13 C. Chandler initial release for CASA 4.1
# 09/23/14 C. Chandler modified to work on CASA 4.2.2, and updated to
#          use Perley-Butler 2013 flux density scale
# 06/08/15 C. Chandler modified for CASA 4.3.1
# 06/08/15 C. Chandler modified for CASA 4.4.0
# 07/13/15 E. Momjian for CASA 4.4.0 split target rflag and final uv
#          plots into two scripts
#          Separate calibrators and targets in final rflag
#          Separate calibrators and targets in statwt
# 07/29/15 E. Momjian force SETJY to create the model column, otherwise
#          calibration will not be correct
# 08/25/15 C. Chandler moved plots to after statwt
# 10/13/15 E. Momjian modified for CASA 4.5.0
#          The use of the real vs. virtual model in setjy is a choice (y/n)
#          at the start of the pipeline
# 02/20/16 E. Momjian modified for CASA 4.5.2
#          Using mstransform based split2 instead of split
# 04/12/16 E. Momjian modified for CASA 4.5.3
# 04/12/16 E. Momjian modified for CASA 4.6.0
#          Mstransform based split2 has been renamed as split.
#          Also using Mstransform based hanningsmooth.
# 10/10/16 E. Momjian modified for CASA 4.7.0
#          Added Amp & phase vs. frequency plots in weblog
# 02/10/17 D. Medlin modified for CASA 4.7.1
#          Added change for CASA version check at start of pipeline
######################################################################

# Change version and date below with each svn commit.  Note changes in the
# .../trunk/doc/CHANGELOG.txt and .../trunk/doc/bugs_features.txt files

import sys
import os
import time

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

time_list = runtiming('startup', 'start')
execfile(pipepath + 'EVLA_pipe_startup.py')
time_list = runtiming('startup', 'end')
pipeline_save()

######################################################################

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
