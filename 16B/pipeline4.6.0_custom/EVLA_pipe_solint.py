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
# DETERMINE SOLINT FOR SCAN-AVERAGE EQUIVALENT

# Split out gain calibrators, drop flagged rows
# NB: in CASA 3.3 gaincal doesn't entirely ignore flagged rows, so for now
# include them in the split below; WILL NEED TO CHANGE TO "keepflags=False"
# FOR CASA 3.4, or use Bryan's method for identifying unflagged data
# in the main ms

logprint ("Starting EVLA_pipe_solint.py", logfileout='logs/solint.log')
time_list=runtiming('solint', 'start')
QA2_solint='Pass'

# Split out calibrators, dropping flagged data; note that it may be
# important *not* to select on field ID, so that fields don't get
# renumbered by split

logprint ("Splitting out calibrators into calibrators.ms", logfileout='logs/solint.log')

syscommand='rm -rf calibrators.ms'
os.system(syscommand)

default('split')
vis=ms_active
outputvis='calibrators.ms'
datacolumn='corrected'
field=''
spw=''
width=int(max(channels))
antenna=''
timebin='0s'
timerange=''
scan=calibrator_scan_select_string
intent=''
array=''
uvrange=''
correlation=''
observation=''
keepflags=False
split()

ms.open('calibrators.ms')
scan_summary = ms.getscansummary()
durations = []
old_spws = []
old_field = ''
for kk in range(len(phase_scan_list)):
    ii = phase_scan_list[kk]

    try:
        end_time = scan_summary[str(ii)]['0']['EndTime']
        begin_time = scan_summary[str(ii)]['0']['BeginTime']
        new_spws = scan_summary[str(ii)]['0']['SpwIds']
        new_field = scan_summary[str(ii)]['0']['FieldId']

        if ((kk > 0) and (phase_scan_list[kk-1] == ii-1) and (set(new_spws) == set(old_spws)) and (new_field == old_field)):
# if contiguous scans, just increase the time on the previous one
            durations[-1] = 86400*(end_time - old_begin_time)
        else:
            durations.append(86400*(end_time - begin_time))
            old_begin_time = begin_time
#        print ii, '  ', durations[-1]
        logprint ("Scan "+str(ii)+" has "+str(durations[-1])+"s on source", logfileout='logs/solint.log')
        ms.reset()
        old_spws = new_spws
        old_field = new_field

    except KeyError:
        logprint ("WARNING: scan "+str(ii)+" is completely flagged and missing from calibrators.ms", logfileout='logs/solint.log')

ms.close()

longsolint = (np.max(durations))*1.01
gain_solint2=str(longsolint)+'s'

# Until we know what the QA criteria are for this script, leave QA2
# set score to "Pass".

logprint ("QA2 score: "+QA2_solint, logfileout='logs/solint.log')
logprint ("Finished EVLA_pipe_solint.py", logfileout='logs/solint.log')
time_list=runtiming('solint', 'end')

pipeline_save()

