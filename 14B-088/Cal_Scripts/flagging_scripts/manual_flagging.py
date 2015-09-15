
import sys

'''
Plot visibility data for each spw to allow for easy manual flags
'''

try:
    vis_name = sys.argv[1]
    field_name = sys.argv[2]
    corrstring = sys.argv[3]
    starting_spw = int(sys.argv[4])
    bp_scan = sys.argv[5]
except IndexError:
    vis_name = raw_input("MS Name? : ")
    field_name = raw_input("Field Name/Number? : ")
    corrstring = raw_input("Corrstring? : ")
    starting_spw = int(raw_input("SPW to start at? : "))
    bp_scan = raw_input("Bandpass scan? : ")

tb.open(vis_name + '/SPECTRAL_WINDOW')
freqs = tb.getcol('REF_FREQUENCY')
nchans = tb.getcol('NUM_CHAN')
tb.close()

spws = range(starting_spw, len(freqs))


for spw_num in spws:
    nchan = nchans[spw_num]

    print "On " + str(spw_num+1) + " of " + str(len(freqs))

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
