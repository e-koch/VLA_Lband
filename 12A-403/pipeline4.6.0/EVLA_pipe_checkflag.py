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

# CHECKING FLAGGING OF BP AND DELAY CALIBRATORS
# use rflag mode of flagdata

logprint ("Starting EVLA_pipe_checkflag.py", logfileout='logs/checkflag.log')
time_list=runtiming('checkflag', 'start')
QA2_checkflag='Pass'

logprint ("Checking RFI flagging of BP and Delay Calibrators", logfileout='logs/checkflag.log')

checkflagfields=''
if (bandpass_field_select_string == delay_field_select_string):
    checkflagfields = bandpass_field_select_string
else:
    checkflagfields = (bandpass_field_select_string+','+
        delay_field_select_string)

# Run only on the bandpass and delay calibrators - testgainscans

default('flagdata')
vis=ms_active
mode='rflag'
field=checkflagfields
correlation='ABS_'+corrstring
scan=testgainscans
ntime='scan'
combinescans=False
datacolumn='corrected'
winsize=3
timedevscale=4.0
freqdevscale=4.0
extendflags=False
action='apply'
display=''
flagbackup=False
savepars=True
flagdata()

#clearstat()

# Until we know what the QA criteria are for this script, leave QA2
# set score to "Pass".

logprint ("QA2 score: "+QA2_checkflag, logfileout='logs/checkflag.log')
logprint ("Finished EVLA_pipe_checkflag.py", logfileout='logs/checkflag.log')
time_list=runtiming('checkflag', 'end')

pipeline_save()


######################################################################
