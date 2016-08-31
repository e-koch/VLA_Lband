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

# CHECKING FLAGGING OF ALL CALIBRATED DATA, INCLUDING TARGET
# use rflag mode of flagdata

logprint ("Starting EVLA_pipe_targetflag.py", logfileout='logs/targetflag.log')
time_list=runtiming('checkflag', 'start')
QA2_targetflag='Pass'

logprint ("Checking RFI flagging of all targets", logfileout='logs/targetflag.log')

# Run on all calibrator scans

default('flagdata')
vis=ms_active
mode='rflag'
field=''
correlation='ABS_'+corrstring
scan=''
intent='*CALIBRATE*'
ntime='scan'
combinescans=False
datacolumn='corrected'
winsize=3
timedevscale=4.0
freqdevscale=4.0
extendflags=False
action='apply'
display=''
flagbackup=True
savepars=True
flagdata()

# Run on all target scans
# Comment out if science target has strong spectral lines
# or set spw to exclude these strong science spectral lines

default('flagdata')
vis=ms_active
mode='rflag'
field=''
correlation='ABS_'+corrstring
scan=''
intent='*TARGET*'
ntime='scan'
combinescans=False
datacolumn='corrected'
winsize=3
timedevscale=4.0
freqdevscale=4.0
extendflags=False
action='apply'
display=''
flagbackup=True
savepars=True
flagdata()



#clearstat()

# Save final version of flags

default('flagmanager')
vis=ms_active
mode='save'
versionname='finalflags'
comment='Final flags saved after calibrations and rflag'
merge='replace'
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/targetflag.log')

# calculate final flag statistics

default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
final_flags = flagdata()

frac_flagged_on_source2 = 1.0-((start_total-final_flags['flagged'])/init_on_source_vis)

logprint ("Final fraction of on-source data flagged = "+str(frac_flagged_on_source2), logfileout='logs/targetflag.log')

if (frac_flagged_on_source2 >= 0.6):
    QA2_targetflag='Fail'

logprint ("QA2 score: "+QA2_targetflag, logfileout='logs/targetflag.log')
logprint ("Finished EVLA_pipe_targetflag.py", logfileout='logs/targetflag.log')
time_list=runtiming('targetflag', 'end')

pipeline_save()


######################################################################
