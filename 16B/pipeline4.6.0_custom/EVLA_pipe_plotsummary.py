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

#MAKING FINAL UV PLOTS ON ALL SOURCES 

logprint ("Starting EVLA_pipe_plotsummary.py", logfileout='logs/plotsummary.log')
time_list=runtiming('plotsummary', 'start')
QA2_plotsummary='Pass'

logprint ("Making final UV plots", logfileout='logs/plotsummary.log')

# Make some plots of the calibrated data

default('plotms')
vis=ms_active
xaxis='time'
yaxis='phase'
ydatacolumn='corrected'
selectdata=True
field=calibrator_field_select_string
correlation=corrstring
averagedata=True
avgchannel=str(max(channels))
avgtime='1e8s'
avgscan=False
transform=False
extendflag=False
iteraxis=''
coloraxis='antenna2'
plotrange=[]
title='Calibrated phase vs. time, all calibrators'
xlabel=''
ylabel=''
showmajorgrid=False
showminorgrid=False
plotfile='all_calibrators_phase_time.png'
overwrite=True
showgui=False
plotms()
#
# see whether we don't need this any more in CASA 4.2.2...
#
#mylogfile = casalog.logfile()
#countmax = 100
#countr = 0
#foundend=False
#while not foundend and countr<countmax:
#    os.system('sleep 10s')
#    f = os.popen('tail -n 10 '+mylogfile)
#    fstrs = f.readlines()
#    f.close()
#    for fstr in fstrs:
#        if fstr.count('End Task: plotms')>0:
#            foundend=True
#            print 'Found end of task plotms in logfile at count '+str(countr)
#    countr+=1

#default('plotms')
#vis=ms_active
#xaxis='time'
#yaxis='amp'
#ydatacolumn='residual'
#selectdata=True
#field=calibrator_field_select_string
#correlation=corrstring
#averagedata=True
#avgchannel=str(max(channels))
#avgtime='1e8s'
#avgscan=False
#transform=False
#extendflag=False
#iteraxis=''
#coloraxis='antenna2'
#plotrange=[]
#title='Corrected-model amp vs. time, all calibrators'
#xlabel=''
#ylabel=''
#showmajorgrid=False
#showminorgrid=False
#plotfile='all_calibrators_resid_amp_time.png'
#overwrite=True
#plotms()
##
#mylogfile = casalog.logfile()
#countmax = 100
#countr = 0
#foundend=False
#while not foundend and countr<countmax:
#    os.system('sleep 10s')
#    f = os.popen('tail -n 10 '+mylogfile)
#    fstrs = f.readlines()
#    f.close()
#    for fstr in fstrs:
#        if fstr.count('End Task: plotms')>0:
#            foundend=True
#            print 'Found end of task plotms in logfile at count '+str(countr)
#    countr+=1

#for ii in range(0,len(calibrator_field_list)):
for ii in field_ids:
    print ii
    default('plotms')
    vis=ms_active
    xaxis='uvwave'
    yaxis='amp'
    ydatacolumn='corrected'
    selectdata=True
#    field=str(calibrator_field_list[ii])
    field=str(field_ids[ii])
    correlation=corrstring
    averagedata=True
    avgchannel=str(max(channels))
    avgtime='1e8s'
    avgscan=False
    transform=False
    extendflag=False
    iteraxis=''
    coloraxis='spw'
    plotrange=[]
    title='Field '+field+', '+field_names[ii]
    xlabel=''
    ylabel=''
    showmajorgrid=False
    showminorgrid=False
    plotfile='field'+field+'_amp_uvdist.png'
    overwrite=True
    showgui=False
    plotms()
#
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
    #Remove plot with no points
    blankstatus = checkblankplot(plotfile, maincasalog)

logprint ("QA2 score: "+QA2_plotsummary, logfileout='logs/plotsummary.log')
logprint ("Finished EVLA_pipe_plotsummary.py", logfileout='logs/plotsummary.log')
time_list=runtiming('plotsummary', 'end')

pipeline_save()


######################################################################
