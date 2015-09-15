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

logprint("Starting EVLA_pipe_flag_uncalspws1b.py", logfileout='logs/flag_uncalspws1b.log')
time_list=runtiming('flag_uncalspws1b', 'start')
QA2_uncalspws1b='Pass'

myscans = scandict

myspw = []
for idd in myscans['DataDescription'].keys():
    ispw = myscans['DataDescription'][idd]['spw']
    if myspw.count(ispw)<1:
        myspw.append(ispw)

calflagresult = getCalFlaggedSoln('delay.k')
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

calflagresult = getCalFlaggedSoln('BPcal.b')
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

flagspw1b = ','.join(["%s" % ii for ii in uniq(flagspwlist)])

if (flagspw1b == ''):
    logprint ("All spws have calibration", logfileout='logs/flag_uncalspws1b.log')
elif (flagspw1b==flagspw1):
    logprint ("No calibration found for spw(s) "+flagspw1b+", already flagged", logfileout='logs/flag_uncalspws1b.log')
else:
    logprint ("No calibration found for spw(s) "+flagspw1b+", flagging these spws in the ms", logfileout='logs/flag_uncalspws1b.log')

    #Now flag with spw=flagspw1b
    default('flagdata')
    vis=ms_active
    mode='manual'
    spw=flagspw1b
    action='apply'
    flagbackup=False
    savepars=True
    async=False
    flagdata()

if (len(uniq(flagspwlist))/float(numSpws2) >= 0.4):
    QA2_uncalspws1b='Fail'
elif (len(uniq(flagspwlist))/float(numSpws2) >= 0.15):
    QA2_uncalspws1b='Partial'

logprint ("QA2 score: "+QA2_uncalspws1b, logfileout='logs/flag_uncalspws1b.log')
logprint ("Finished EVLA_pipe_flag_uncalspws1b.py", logfileout='logs/flag_uncalspws1b.log')
time_list=runtiming('flag_uncalspws1b', 'end')

pipeline_save()
