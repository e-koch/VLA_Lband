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

# CHECK THE FLAGS BEFORE FINAL CALIBRATION
# APPLY ALL CALIBRATIONS AND CHECK CALIBRATED DATA
# CHECK THE FLAGS AFTER FINAL CALIBRATION


logprint ("Starting EVLA_pipe_applycals.py", logfileout='logs/applycals.log')
time_list=runtiming('applycals', 'start')
QA2_applycals='Pass'

default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
myinitialflags = flagdata()
logprint ("Finished flags summary before final applycal", logfileout='logs/applycals.log')


FinalGainTables=copy.copy(priorcals)
FinalGainTables.append('finaldelay.k')
FinalGainTables.append('finalBPcal.b')
FinalGainTables.append('averagephasegain.g')
FinalGainTables.append('finalampgaincal.g')
FinalGainTables.append('finalphasegaincal.g')

ntables=len(FinalGainTables)

default('applycal')
vis=ms_active
field=''
spw=''
intent=''
selectdata=False
gaintable=FinalGainTables
gainfield=['']
interp=['']
spwmap=[]
parang=False
calwt=[False]*ntables
applymode='calflagstrict'
flagbackup=True
applycal()


default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
myinitialflags = flagdata()
logprint ("Finished flags summary after final applycal", logfileout='logs/applycals.log')


# Check calibrated data for more flagging

# Plot "corrected" datacolumn for calibrators in plotms
# Plot "corrected" datacolumn for targets
# Amp vs. freq, amp vs. time

# Until we understand better the failure modes of this task, leave
# QA2 score set to "Pass".

logprint ("QA2 score: "+QA2_applycals, logfileout='logs/applycals.log')
logprint ("Finished EVLA_pipe_applycals.py", logfileout='logs/applycals.log')
time_list=runtiming('applycals', 'end')

pipeline_save()
