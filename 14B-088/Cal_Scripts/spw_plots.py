
'''
Plot UVdist per field per SPW
Plot BP amp and phase per antenna per SPW

Requires that the pipeline namespace be populated.
'''

import os
import numpy as np
import sys

# Repopulate namespace
execfile("/home/ekoch/canfar_scripts/EVLA_pipeline1.3.0/EVLA_pipe_restore.py")

print 'Now plotting SPW plots...'
print ms_active
print pipepath

# Specify the type of plots which are created
# In the command line, specify 'T' or 'F'
uv_plots = sys.argv[4]
if uv_plots == 'T':
    uv_plots = True
else:
    uv_plots = False
ampphase_time_plots = sys.argv[5]
if ampphase_time_plots == 'T':
    ampphase_time_plots = True
else:
    ampphase_time_plots = False

bpcal_plots = sys.argv[6]
if bpcal_plots == 'T':
    bpcal_plots = True
else:
    bpcal_plots = False

# UV plots per SPW
if uv_plots:
    print "Creating UV plots per field per SPW, color by antenna"
    for ii in field_ids:
        for jj in field_spws[ii]:
            print ii, jj
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
            coloraxis = 'ant2'
            plotrange = []
            title = 'Field ' + field + ', ' + field_names[ii] + " SPW " + spw
            xlabel = ''
            ylabel = ''
            showmajorgrid = False
            showminorgrid = False
            plotfile = 'field' + field + '_SPW_' + spw + '_amp_uvdist.png'
            overwrite = True
            showgui = False
            async = False
            plotms()

if ampphase_time_plots:
    print "Creating Amp vs Time per field, color by antenna"
    for ii in field_ids:
        print ii
        default('plotms')
        vis = ms_active
        xaxis = 'time'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = str(ii)
        spw = ''
        correlation = corrstring
        averagedata = True
        avgchannel = str(max(channels))
        avgtime = '1e8s'
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'ant2'
        plotrange = []
        title = 'Field ' + field + ', ' + field_names[ii]
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotfile = 'field' + field + '_amp_time.png'
        overwrite = True
        showgui = False
        async = False
        plotms()

    print "Creating Phase vs Time per field, color by antenna"

    for ii in field_ids:
        print ii
        default('plotms')
        vis = ms_active
        xaxis = 'time'
        yaxis = 'phase'
        ydatacolumn = 'corrected'
        selectdata = True
        field = str(ii)
        spw = ''
        correlation = corrstring
        averagedata = True
        avgchannel = str(max(channels))
        avgtime = '1e8s'
        avgscan = False
        transform = False
        extendflag = False
        iteraxis = ''
        coloraxis = 'ant2'
        plotrange = []
        title = 'Field ' + field + ', ' + field_names[ii]
        xlabel = ''
        ylabel = ''
        showmajorgrid = False
        showminorgrid = False
        plotfile = 'field' + field + '_phase_time.png'
        overwrite = True
        showgui = False
        async = False
        plotms()

# Final BP plots per SPW
if bpcal_plots:

    # Somewhere in the cal tables (though not accessible)

    numSpws = len(field_spws[0])

    nplots = int(numSpws)

    # if ((numSpws % 3) > 0):
    #     nplots = nplots + 1

    # To plot the cal tables, we need to move to the same directory as the MS
    ms_path = "/".join(ms_active.split('/')[:-1])
    ms_name = ms_active.split('/')[-1]

    proc_path = os.getcwd()

    os.chdir(ms_path)

    # Need to append the correct name to the MS
    tb.open('finalBPcal.b/', nomodify=False,
            lockoptions={'option': 'user'})

    tb.lock()
    tb.removekeyword('MSName')
    tb.putkeyword('MSName', ms_name)
    tb.unlock()

    dataVarCol = tb.getvarcol('CPARAM')
    flagVarCol = tb.getvarcol('FLAG')
    tb.close()

    tb.open('finalBPcal.b/')
    print "MSName in caltable changed to:" + str(tb.getkeyword('MSName'))
    tb.close()

    rowlist = dataVarCol.keys()
    nrows = len(rowlist)
    maxmaxamp = 0.0
    maxmaxphase = 0.0
    for rrow in rowlist:
        dataArr = dataVarCol[rrow]
        flagArr = flagVarCol[rrow]
        amps = np.abs(dataArr)
        phases = np.arctan2(np.imag(dataArr), np.real(dataArr))
        good = np.logical_not(flagArr)
        tmparr = amps[good]
        if (len(tmparr) > 0):
            maxamp = np.max(amps[good])
            if (maxamp > maxmaxamp):
                maxmaxamp = maxamp
        tmparr = np.abs(phases[good])
        if (len(tmparr) > 0):
            maxphase = np.max(np.abs(phases[good])) * 180. / np.pi
            if (maxphase > maxmaxphase):
                maxmaxphase = maxphase
    ampplotmax = maxmaxamp
    phaseplotmax = maxmaxphase

    for ii in range(1, nplots+1):
        filename = proc_path+'/finalBPcal_amp' + str(ii) + '.png'
        syscommand = 'rm -rf ' + filename
        os.system(syscommand)

        spwPlot = str(ii) # + '~' + str(ii * 3 + 2)

        default('plotcal')
        caltable = 'finalBPcal.b'
        xaxis = 'freq'
        yaxis = 'amp'
        poln = ''
        field = ''
        antenna = ''
        spw = spwPlot
        timerange = ''
        subplot = 111
        overplot = False
        clearpanel = 'Auto'
        iteration = 'spw'
        plotrange = [0, 0, 0, ampplotmax]
        showflags = False
        plotsymbol = 'o'
        plotcolor = 'blue'
        markersize = 5.0
        fontsize = 10.0
        showgui = False
        figfile = filename
        async = False
        plotcal()

    for ii in range(1, nplots+1):
        filename = proc_path+'/finalBPcal_phase' + str(ii) + '.png'
        syscommand = 'rm -rf ' + filename
        os.system(syscommand)

        spwPlot = str(ii) # + '~' + str(ii * 3 + 2)

        default('plotcal')
        caltable = 'finalBPcal.b'
        xaxis = 'freq'
        yaxis = 'phase'
        poln = ''
        field = ''
        antenna = ''
        spw = spwPlot
        timerange = ''
        subplot = 111
        overplot = False
        clearpanel = 'Auto'
        iteration = 'spw'
        plotrange = [0, 0, -phaseplotmax, phaseplotmax]
        showflags = False
        plotsymbol = 'o'
        plotcolor = 'blue'
        markersize = 5.0
        fontsize = 10.0
        showgui = False
        figfile = filename
        async = False
        plotcal()

    os.chdir(proc_path)
