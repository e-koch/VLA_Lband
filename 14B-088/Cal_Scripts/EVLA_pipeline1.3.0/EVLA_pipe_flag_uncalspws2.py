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

# Determine uncalibrated spws in the MS and flag them
# Steve Myers notes that this should be handled in the long term by applycal
# with a new option applymode='strict'

logprint("Starting EVLA_pipe_flag_uncalspws2.py", logfileout='logs/flag_uncalspws2.log')
time_list=runtiming('flag_uncalspws2', 'start')
QA2_uncalspws2='Pass'

myscans = scandict

myspw = []
for idd in myscans['DataDescription'].keys():
    ispw = myscans['DataDescription'][idd]['spw']
    if myspw.count(ispw)<1:
        myspw.append(ispw)

calflagresult = getCalFlaggedSoln('finaldelay.k')
goodspw = []
for ispw in calflagresult['spw'].keys():
    tot = 0.0
    flagd = 0.0
    for ipol in calflagresult['spw'][ispw].keys():
        tot += calflagresult['spw'][ispw][ipol]['total']
        flagd += calflagresult['spw'][ispw][ipol]['flagged']
    if tot>0:
        fract = flagd/tot
        if fract<1.0:
            goodspw.append(ispw)

flagspwlist = []
flagspw = ''
for ispw in myspw:
    if goodspw.count(ispw)<1:
        flagspwlist.append(ispw)
        if flagspw=='':
            flagspw = str(ispw)
        else:
            flagspw += ','+str(ispw)

calflagresult = getCalFlaggedSoln('finalBPcal.b')
goodspw = []
for ispw in calflagresult['spw'].keys():
    tot = 0.0
    flagd = 0.0
    for ipol in calflagresult['spw'][ispw].keys():
        tot += calflagresult['spw'][ispw][ipol]['total']
        flagd += calflagresult['spw'][ispw][ipol]['flagged']
    if tot>0:
        fract = flagd/tot
        if fract<1.0:
            goodspw.append(ispw)

for ispw in myspw:
    if goodspw.count(ispw)<1:
        flagspwlist.append(ispw)
        if flagspw=='':
            flagspw = str(ispw)
        else:
            flagspw += ','+str(ispw)

flagspw2 = ','.join(["%s" % ii for ii in uniq(flagspwlist)])


if (flagspw2 == ''):
    logprint ("All spws have calibration", logfileout='logs/flag_uncalspws2.log')
elif (flagspw2==flagspw1b):
    logprint ("No calibration found for spw(s) "+flagspw2+", already flagged", logfileout='logs/flag_uncalspws2.log')
else:
    logprint ("No calibration found for spw(s) "+flagspw2+", flagging these spws in the ms", logfileout='logs/flag_uncalspws2.log')

    #Now flag with spw=flagspw2
    default('flagdata')
    vis=ms_active
    mode='manual'
    spw=flagspw2
    action='apply'
    flagbackup=False
    savepars=True
    async=False
    flagdata()

if (len(uniq(flagspwlist))/float(numSpws2) >= 0.4):
    QA2_uncalspws2='Fail'
elif (len(uniq(flagspwlist))/float(numSpws2) >= 0.15):
    QA2_uncalspws2='Partial'

logprint ("QA2 score: "+QA2_uncalspws2, logfileout='logs/flag_uncalspws2.log')
logprint ("Finished EVLA_pipe_flag_uncalspws2.py", logfileout='logs/flag_uncalspws2.log')
time_list=runtiming('flag_uncalspws2', 'end')

pipeline_save()
