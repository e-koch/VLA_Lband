
import sys
from glob import glob
import os
from datetime import datetime

from tasks import plotms

'''
Plot visibility data for each spw to allow for easy manual flags
'''

try:
    vis_name = sys.argv[1]
    field_names = sys.argv[2]
    corrstring = sys.argv[3]
    starting_spw = int(sys.argv[4])
    bp_scan = sys.argv[5]
except IndexError:
    vis_name = raw_input("MS Name? : ")
    field_names = raw_input("Field Name/Number(s)? : ")
    corrstring = raw_input("Corrstring? : ")
    starting_spw = int(raw_input("SPW to start at? : "))
    bp_scan = raw_input("Bandpass scan? : ")

tb.open(vis_name + '/SPECTRAL_WINDOW')
freqs = tb.getcol('REF_FREQUENCY')
nchans = tb.getcol('NUM_CHAN')
tb.close()

spws = range(starting_spw, len(freqs))

fields = field_names.split(",")

for n, field_name in enumerate(fields):

    print("On {0}. {1} out of {2} fields.".format(field_name, n + 1,
                                                  len(fields)))

    for spw_num in spws:
        nchan = nchans[spw_num]

        print "On SPW {0} of {1}".format(str(spw_num + 1), str(len(freqs)))

        default('plotms')
        vis = vis_name
        xaxis = 'time'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = field_name
        spw = str(spw_num)
        scan = bp_scan
        correlation = corrstring
        averagedata = False
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'antenna2'
        plotrange = []
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotms()

        raw_input("Continue?")

        default('plotms')
        vis = vis_name
        xaxis = 'channel'
        yaxis = 'phase'
        ydatacolumn = 'corrected'
        selectdata = True
        field = field_name
        spw = str(spw_num)
        correlation = corrstring
        averagedata = True
        avgtime = '1e8s'
        avgscan = True
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'antenna2'
        plotrange = []
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotms()

        raw_input("Continue?")

        default('plotms')
        vis = vis_name
        xaxis = 'channel'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = field_name
        spw = str(spw_num)
        correlation = corrstring
        averagedata = True
        avgtime = '1e8s'
        avgscan = True
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'antenna2'
        plotrange = []
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotms()

        raw_input("Continue?")

        default('plotms')
        vis = vis_name
        xaxis = 'time'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = field_name
        spw = str(spw_num)
        correlation = corrstring
        averagedata = True
        avgchannel = str(nchan)
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = 'scan'
        coloraxis = 'antenna2'
        plotrange = []
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotms()

        raw_input("Continue?")

        default('plotms')
        vis = vis_name
        xaxis = 'time'
        yaxis = 'phase'
        ydatacolumn = 'corrected'
        selectdata = True
        field = field_name
        spw = str(spw_num)
        correlation = corrstring
        averagedata = True
        avgchannel = str(nchan)
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = 'scan'
        coloraxis = 'antenna2'
        plotrange = []
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotms()

        raw_input("Continue?")

        default('plotms')
        vis = vis_name
        xaxis = 'uvwave'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = field_name
        spw = str(spw_num)
        correlation = corrstring
        averagedata = True
        avgchannel = str(nchan)
        avgtime = '1e8s'
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'antenna2'
        plotrange = []
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotms()

        raw_input("Continue?")

# Get the existing flag version names.
flag_folder = "{}.flagversions".format(vis_name)
tstamp = datetime.now().strftime("%Y%m%d-%H%M%S")
if not os.path.exists(flag_folder):
    print("No flag versions exist. Using default flag name.")
    versionname = "manual_flagging_1_{}".format(tstamp)
else:
    flag_versions = glob(os.path.join(flag_folder, "flag.manual_flagging_*"))
    if len(flag_versions) == 0:
        versionname = "manual_flagging_1_{}".format(tstamp)
    else:
        num = len(flag_versions) + 1
        versionname = "manual_flagging_{0}_{1}".format(num, tstamp)

# Save this new version of the flags
flagmanager(vis=vis_name, mode='save',
            versionname=versionname)
