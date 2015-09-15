
'''
Produces the plots made with plotms during the pipeline without recomputing.
'''

pipepath = '/lustre/aoc/observers/nm-7669/canfar_scripts/EVLA_pipeline1.3.0/'

execfile(pipepath+"EVLA_pipe_restore.py")

print("Plotting all calibrators phase vs. time")

default('plotms')
vis = ms_active
xaxis = 'time'
yaxis = 'phase'
ydatacolumn = 'corrected'
selectdata = True
field = calibrator_field_select_string
correlation = corrstring
averagedata = True
avgchannel = str(max(channels))
avgtime = '1e8s'
avgscan = False
transform = False
extendflag = False
iteraxis = ''
coloraxis = 'antenna2'
plotrange = []
title = 'Calibrated phase vs. time, all calibrators'
xlabel = ''
ylabel = ''
showmajorgrid = False
showminorgrid = False
plotfile = 'all_calibrators_phase_time.png'
overwrite = True
showgui = False
async = False
plotms()

print("Plotting individual SPWs by calibrator, phase vs time and amp vs time")

for ii in calibrator_field_select_string.split(","):
    for jj in field_spws[0]:

        print("Calibrator: %s SPW: %s" % (ii, jj))

        default('plotms')
        vis = ms_active
        xaxis = 'time'
        yaxis = 'phase'
        ydatacolumn = 'corrected'
        selectdata = True
        field = str(ii)
        spw = str(jj)
        correlation = corrstring
        averagedata = True
        avgchannel = str(max(channels))
        avgtime = '1e8s'
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'antenna2'
        plotrange = []
        title = 'Calibrated phase vs. time, calibrator' + str(ii)
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotfile = 'calibrator_' + \
            str(ii) + '_spw_' + str(jj) + '_phase_time.png'
        overwrite = True
        showgui = False
        async = False
        plotms()

        default('plotms')
        vis = ms_active
        xaxis = 'time'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = str(ii)
        spw = str(jj)
        correlation = corrstring
        averagedata = True
        avgchannel = str(max(channels))
        avgtime = '1e8s'
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'antenna2'
        plotrange = []
        title = 'Calibrated amp vs. time, calibrator' + str(ii)
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotfile = 'calibrator_' + \
            str(ii) + '_spw_' + str(jj) + '_amp_time.png'
        overwrite = True
        showgui = False
        async = False
        plotms()


print("Plotting UVWave vs Amplitude.")

for ii in field_ids:
    print("On Field: " + str(ii))
    default('plotms')
    vis = ms_active
    xaxis = 'uvwave'
    yaxis = 'amp'
    ydatacolumn = 'corrected'
    selectdata = True
#    field=str(calibrator_field_list[ii])
    field = str(field_ids[ii])
    correlation = corrstring
    averagedata = True
    avgchannel = str(max(channels))
    avgtime = '1e8s'
    avgscan = False
    transform = False
    extendflag = False
    iteraxis = ''
    coloraxis = 'spw'
    plotrange = []
    title = 'Field ' + field + ', ' + field_names[ii]
    xlabel = ''
    ylabel = ''
    showmajorgrid = False
    showminorgrid = False
    plotfile = 'field' + field + '_amp_uvdist.png'
    overwrite = True
    showgui = False
    async = False
    plotms()

print("Plotting UVWave vs Amp. per SPW.")

for ii in field_ids:
    for jj in field_spws[ii]:
        print("Field: %s SPW: %s" % (str(ii), str(jj)))

        default('plotms')
        vis = ms_active
        xaxis = 'uvwave'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = str(ii)
        spw = str(jj)
        correlation = corrstring
        averagedata = True
        avgchannel = str(max(channels))
        avgtime = '1e8s'
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'antenna2'
        plotrange = []
        title = 'Field ' + field + ', ' + field_names[ii] + " SPW " + spw
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotfile = 'field' + field + '_SPW_' + spw + '_amp_uvwave.png'
        overwrite = True
        showgui = False
        async = False
        plotms()
