
import numpy as np
import re
import os
import sys

if len(sys.argv) > 3:
    ms_active = sys.argv[-3]
    field_str = sys.argv[-2]
    spw_str = sys.argv[-1]
else:
    ms_active = raw_input("MS? : ")
    field_str = raw_input("Field? : ")
    spw_str = raw_input("SPW? : ")


tb.open(ms_active+"/FIELD")
names = tb.getcol('NAME')
matches = [string for string in names if re.match(field_str, string)]
posn_matches = \
    [i for i, string in enumerate(names) if re.match(field_str, string)]
numFields = tb.nrows()
tb.close()

if len(matches) == 0:
    raise TypeError("No matches found for the given field string")

tb.open(ms_active)
scanNums = sorted(np.unique(tb.getcol('SCAN_NUMBER')))
field_scans = []
for ii in range(numFields):
    subtable = tb.query('FIELD_ID==%s'%ii)
    field_scans.append(list(np.unique(subtable.getcol('SCAN_NUMBER'))))
tb.close()

field_scans = [scans for i, scans in enumerate(field_scans) if i in posn_matches]

scan_dir = ms_active.rstrip(".ms")+"_scan_plots"

try:
    os.mkdir(scan_dir)
except OSError:
    pass

for ii in range(len(field_scans)):
    print("On field "+matches[ii])
    for jj in field_scans[ii]:

        print("On scan "+str(jj))

        default('plotms')
        vis = ms_active
        xaxis = 'time'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = matches[ii]
        scan = str(jj)
        spw = spw_str
        correlation = "RR,LL"
        averagedata = True
        avgbaseline = True
        transform = False
        extendflag = False
        plotrange = []
        title = 'Amp vs Time: Field '+matches[ii]+' Scan'+str(jj)
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotfile = scan_dir+'/field_'+matches[ii]+'_scan_'+str(jj)+'.png'
        overwrite = True
        showgui = False
        async = False
        plotms()
