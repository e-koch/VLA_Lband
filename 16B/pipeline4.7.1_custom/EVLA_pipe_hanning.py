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

# HANNING SMOOTH (OPTIONAL, IMPORTANT IF THERE IS NARROWBAND RFI)
# NB: Steve points out that really you want to do Hanning smoothing
# only after a delay calibration, to avoid decorrelation
# NB: myHanning is set to "n" by default, but this should be enabled
# for the user to choose during reprocessing
# 04/13/16 E. Momjian: hanning smoothing is now mstransform based
#                      and does not allow overwriting the input MS.
#                      Using B.K.'s code to replace the input MS with
#                      the task's output MS.
#                      Also, copying *.xml to the output MS for 
#                      flagcmd not to fail as it needs Flags.xml

import shutil
import glob
import os

logprint ("Starting EVLA_pipe_hanning.py", logfileout='logs/hanning.log')
time_list=runtiming('hanning', 'start')
QA2_hanning='Pass'

if myHanning.lower() == "y":
#    if (os.path.exists(mshsmooth) == False):
        logprint ("Hanning smoothing the data", logfileout='logs/hanning.log')

        default('hanningsmooth')
        vis=msname
        datacolumn='data'
#        outputvis=mshsmooth
        outputvis='temphanning.ms'
        hanningsmooth()
        myHanning="n"

        logprint ('Copying xml files to the output ms')
        for file in glob.glob(msname+'/*.xml'):
                shutil.copy2(file , 'temphanning.ms/')
        logprint ('Removing original VIS '+msname, logfileout='logs/hanning.log')
        shutil.rmtree(msname)
        logprint('Renaming temphanning.ms to '+msname, logfileout='logs/hanning.log')
        os.rename('temphanning.ms', msname)
        logprint ("Hanning smoothing finished, myHanning parameter reset to 'n' to avoid further smoothing on restarts", logfileout='logs/hanning.log')
#    else:
#        logprint ("Hanning smoothed ms set already exists, will use existing ms", logfileout='logs/hanning.log')
else:
    logprint ("NOT Hanning smoothing the data", logfileout='logs/hanning.log')


# Until we know better the possible failures modes of this script,
# leave set QA2 score set to "Pass".

logprint ("Finished EVLA_pipe_hanning.py", logfileout='logs/hanning.log')
logprint ("QA2 score: "+QA2_hanning, logfileout='logs/hanning.log')

time_list=runtiming('hanning', 'end')

pipeline_save()
