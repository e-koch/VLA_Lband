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

logprint("Starting EVLA_pipe_targetflag.py", logfileout='logs/targetflag.log')
time_list = runtiming('checkflag', 'start')
QA2_targetflag = 'Pass'

logprint("Checking RFI flagging of all targets",
         logfileout='logs/targetflag.log')

# Run on all targets/scans

default('flagdata')
vis = ms_active
mode = 'rflag'
field = ''
correlation = 'ABS_' + corrstring
scan = ''
ntime = 'scan'
combinescans = False
datacolumn = 'corrected'
winsize = 3
timedevscale = 4.0
freqdevscale = 4.0
extendflags = False
action = 'apply'
display = ''
flagbackup = True
savepars = True
async = False
flagdata()


# clearstat()

# Save final version of flags

default('flagmanager')
vis = ms_active
mode = 'save'
versionname = 'finalflags'
comment = 'Final flags saved after calibrations and rflag'
merge = 'replace'
async = False
flagmanager()
logprint("Flag column saved to " + versionname,
         logfileout='logs/targetflag.log')

# Make some plots of the calibrated data

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
title = 'Calibrated phase vs. time, calibrator' + str(ii)
xlabel = ''
ylabel = ''
showmajorgrid = False
showminorgrid = False
plotfile = 'all_calibrators_phase_time.png'
overwrite = True
showgui = False
async = False
plotms()

for ii in calibrator_field_select_string.split(","):
    for jj in field_spws[0]:

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

#
# see whether we don't need this any more in CASA 4.2.2...
#
#mylogfile = casalog.logfile()
#countmax = 100
#countr = 0
# foundend=False
# while not foundend and countr<countmax:
#    os.system('sleep 10s')
#    f = os.popen('tail -n 10 '+mylogfile)
#    fstrs = f.readlines()
#    f.close()
#    for fstr in fstrs:
#        if fstr.count('End Task: plotms')>0:
#            foundend=True
#            print 'Found end of task plotms in logfile at count '+str(countr)
#    countr+=1

# default('plotms')
# vis=ms_active
# xaxis='time'
# yaxis='amp'
# ydatacolumn='residual'
# selectdata=True
# field=calibrator_field_select_string
# correlation=corrstring
# averagedata=True
# avgchannel=str(max(channels))
# avgtime='1e8s'
# avgscan=False
# transform=False
# extendflag=False
# iteraxis=''
# coloraxis='antenna2'
# plotrange=[]
#title='Corrected-model amp vs. time, all calibrators'
# xlabel=''
# ylabel=''
# showmajorgrid=False
# showminorgrid=False
# plotfile='all_calibrators_resid_amp_time.png'
# overwrite=True
# async=False
# plotms()
##
#mylogfile = casalog.logfile()
#countmax = 100
#countr = 0
# foundend=False
# while not foundend and countr<countmax:
#    os.system('sleep 10s')
#    f = os.popen('tail -n 10 '+mylogfile)
#    fstrs = f.readlines()
#    f.close()
#    for fstr in fstrs:
#        if fstr.count('End Task: plotms')>0:
#            foundend=True
#            print 'Found end of task plotms in logfile at count '+str(countr)
#    countr+=1

# for ii in range(0, len(calibrator_field_list)):
for ii in field_ids:
    print ii
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
        coloraxis = 'antenna2'
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


# see if we can omit this now in CASA 4.2.2...
#
#    mylogfile = casalog.logfile()
#    countmax = 100
#    countr = 0
#    foundend=False
#    while not foundend and countr<countmax:
#        os.system('sleep 10s')
#        f = os.popen('tail -n 10 '+mylogfile)
#        fstrs = f.readlines()
#        f.close()
#        for fstr in fstrs:
#            if fstr.count('End Task: plotms')>0:
#                foundend=True
#                print 'Found end of task plotms in logfile at count '+str(countr)
#        countr+=1
#
    # Remove plot with no points
#    blankstatus = checkblankplot(plotfile, maincasalog)

# calculate final flag statistics

default('flagdata')
vis = ms_active
mode = 'summary'
spwchan = True
spwcorr = True
basecnt = True
action = 'calculate'
savepars = False
async = False
final_flags = flagdata()

frac_flagged_on_source2 = 1.0 - \
    ((start_total - final_flags['flagged']) / init_on_source_vis)

logprint("Final fraction of on-source data flagged = " +
         str(frac_flagged_on_source2), logfileout='logs/targetflag.log')

if (frac_flagged_on_source2 >= 0.6):
    QA2_targetflag = 'Fail'

logprint("QA2 score: " + QA2_targetflag, logfileout='logs/targetflag.log')
logprint("Finished EVLA_pipe_targetflag.py", logfileout='logs/targetflag.log')
time_list = runtiming('targetflag', 'end')


pipeline_save()


######################################################################
