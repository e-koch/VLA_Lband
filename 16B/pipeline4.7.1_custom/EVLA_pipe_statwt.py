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

# CALCULATE DATA WEIGHTS BASED ON ST. DEV. WITHIN EACH SPW
# use statwt

logprint ("Starting EVLA_pipe_statwt.py", logfileout='logs/statwt.log')
time_list=runtiming('checkflag', 'start')
QA2_statwt='Pass'

logprint ("Calculate data weights per spw using statwt", logfileout='logs/statwt.log')

# Run on all calibrators

default(statwt)
vis=ms_active
dorms=False
fitspw=''
fitcorr=''
combine=''
minsamp=2
field=''
spw=''
intent='*CALIBRATE*'
datacolumn='corrected'
statwt()

# Run on all targets
# set spw to exclude strong science spectral lines

default(statwt)
vis=ms_active
dorms=False
fitspw=''
fitcorr=''
combine=''
minsamp=2
field=''
spw=''
intent='*TARGET*'
datacolumn='corrected'
statwt()

# Until we understand better the failure modes of this task, leave QA2
# score set to "Pass".

logprint ("QA2 score: "+QA2_statwt, logfileout='logs/statwt.log')
logprint ("Finished EVLA_pipe_statwt.py", logfileout='logs/statwt.log')
time_list=runtiming('targetflag', 'end')

pipeline_save()


######################################################################
