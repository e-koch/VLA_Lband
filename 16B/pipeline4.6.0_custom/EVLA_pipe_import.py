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

# IMPORT THE DATA TO CASA

# NB: apply shadow and zero flags, but apply online flags later when
# tbuff has been determined from the integration time
# Version 20121108 no flags generated here, do in EVLA_pipe_flagall.py

logprint ("Starting EVLA_pipe_import.py", logfileout='logs/import.log')
time_list=runtiming('import', 'start')
QA2_import='Pass'

if (os.path.exists(msname) == False):
    logprint ("Creating measurementset", logfileout='logs/import.log')

    default('importevla')
    asdm=SDM_name
    vis=msname
    ocorr_mode='co'
    compression=False
    asis=''
    scans=''
    verbose=True
    overwrite=False
    online=False
    flagzero=False
    flagpol=False
    shadow=False
    tolerance=0.0
    addantenna=''
    applyflags=False
    savecmds=False
    flagbackup=False
    importevla()

    logprint ("Measurement set "+msname+" created", logfileout='logs/import.log')

else:
    logprint ("Measurement set already exists, will use "+msname, logfileout='logs/import.log')

# Until we understand better the possible failure modes to look for
# in this script, leave QA2 set to "Pass".

logprint ("Finished EVLA_pipe_import.py", logfileout='logs/import.log')
logprint ("QA2 score: "+QA2_import, logfileout='logs/import.log')
time_list=runtiming('import', 'end')

pipeline_save()

