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

# SEMI-FINAL DELAY AND BANDPASS CALIBRATIONS
# (semi-final because we have not yet determined the spectral index 
# of the bandpass calibrator)

# Find reference antenna again after all that flagging

logprint ("Starting EVLA_pipe_semiFinalBPdcals1.py", logfileout='logs/semiFinalBPdcals1.log')
time_list=runtiming('semiFinalBPdcals1', 'start')
QA2_semiFinalBPdcals1='Pass'

logprint ("Finding a reference antenna for semi-final delay and BP calibrations", logfileout='logs/semiFinalBPdcals1.log')


refantspw=''
refantfield=calibrator_field_select_string

findrefant=RefAntHeuristics(vis=ms_active,field=refantfield,geometry=True,flagging=True)
RefAntOutput=findrefant.calculate()
refAnt=str(RefAntOutput[0])

logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/semiFinalBPdcals1.log')

# Initial phase solutions on delay calibrator

syscommand='rm -rf semiFinaldelayinitialgain.g'
os.system(syscommand)

if (cal3C84_d == True):
    default('gaincal')
    vis=ms_active
    caltable='semiFinaldelayinitialgain.g'
    field=delay_field_select_string
    spw=tst_delay_spw
    intent=''
    selectdata=True
    uvrange=uvrange3C84
    scan=delay_scan_select_string
    solint='int'
    combine='scan'
    preavg=-1.0
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=3.0
    solnorm=False
    gaintype='G'
    smodel=[]
    calmode='p'
    append=False
    docallib=False
    gaintable=priorcals
    gainfield=['']
    interp=['']
    spwmap=[]
    parang=False
    async=False
    gaincal()
else:
    default('gaincal')
    vis=ms_active
    caltable='semiFinaldelayinitialgain.g'
    field=delay_field_select_string
    spw=tst_delay_spw
    intent=''
    selectdata=True
    uvrange=''
    scan=delay_scan_select_string
    solint='int'
    combine='scan'
    preavg=-1.0
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=3.0
    solnorm=False
    gaintype='G'
    smodel=[]
    calmode='p'
    append=False
    docallib=False
    gaintable=priorcals
    gainfield=['']
    interp=['']
    spwmap=[]
    parang=False
    async=False
    gaincal()

syscommand='rm -rf delay.k'
os.system(syscommand)

flaggedSolnResult=semiFinaldelays(ms_active,'delay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/semiFinalBPdcals1.log')
logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/semiFinalBPdcals1.log')

if (flaggedSolnResult['all']['total'] > 0):
    fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
else:
    fracFlaggedSolns=1.0

if (fracFlaggedSolns > critfrac):
    logprint ("Not enough good solutions, trying a different reference antenna", logfileout='logs/semiFinalBPdcals1.log')
    refAnt=str(RefAntOutput[1])
    logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/semiFinalBPdcals1.log')
    flaggedSolnResult=semiFinaldelays(ms_active,'delay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
    logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/semiFinalBPdcals1.log')
    logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/semiFinalBPdcals1.log')

    if (flaggedSolnResult['all']['total'] > 0):
        fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
    else:
        fracFlaggedSolns=1.0

    if (fracFlaggedSolns > critfrac):
        logprint ("Not enough good solutions, trying a different reference antenna", logfileout='logs/semiFinalBPdcals1.log')
        refAnt=str(RefAntOutput[2])
        logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/semiFinalBPdcals1.log')

        flaggedSolnResult=semiFinaldelays(ms_active,'delay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
        logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/semiFinalBPdcals1.log')
        logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/semiFinalBPdcals1.log')

        if (flaggedSolnResult['all']['total'] > 0):
            fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
        else:
            fracFlaggedSolns=1.0

        if (fracFlaggedSolns > critfrac):
            logprint ("Not enough good solutions, trying a different reference antenna", logfileout='logs/semiFinalBPdcals1.log')
            refAnt=str(RefAntOutput[3])
            logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/semiFinalBPdcals1.log')
            flaggedSolnResult=semiFinaldelays(ms_active,'delay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
            logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/semiFinalBPdcals1.log')
            logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/semiFinalBPdcals1.log')

            if (flaggedSolnResult['all']['total'] > 0):
                fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
            else:
                fracFlaggedSolns=1.0

            if (fracFlaggedSolns > critfrac):
                logprint ("Warning, tried several reference antennas, there might be something wrong with your data", logfileout='logs/semiFinalBPdcals1.log')


logprint ("Delay calibration complete", logfileout='logs/semiFinalBPdcals1.log')

logprint ("Plotting delays", logfileout='logs/semiFinalBPdcals1.log')

nplots=int(numAntenna/3)

if ((numAntenna%3)>0):
    nplots = nplots + 1

for ii in range(nplots):
    filename='delay'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='delay.k'
    xaxis='freq'
    yaxis='delay'
    poln=''
    field=''
    antenna=antPlot
    spw=''
    timerange=''
    subplot=311
    overplot=False
    clearpanel='Auto'
    iteration='antenna'
    plotrange=[]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()

# Do initial gaincal on BP calibrator then semi-final BP calibration

syscommand='rm -rf BPinitialgain.g'
os.system(syscommand)

GainTables=copy.copy(priorcals)
GainTables.append('delay.k')

if (cal3C84_bp == True):
    default('gaincal')
    vis=ms_active
    caltable='BPinitialgain.g'
    field=''
    spw=tst_bpass_spw
    selectdata=True
    uvrange=uvrange3C84
    scan=bandpass_scan_select_string
    solint=gain_solint1
    combine='scan'
    preavg=-1.0
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=3.0
    solnorm=False
    gaintype='G'
    smodel=[]
    calmode='p'
    append=False
    docallib=False
    gaintable=GainTables
    gainfield=['']
    interp=['']
    spwmap=[]
    parang=False
    async=False
    gaincal()
else:
    default('gaincal')
    vis=ms_active
    caltable='BPinitialgain.g'
    field=''
    spw=tst_bpass_spw
    selectdata=True
    uvrange=''
    scan=bandpass_scan_select_string
    solint=gain_solint1
    combine='scan'
    preavg=-1.0
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=3.0
    solnorm=False
    gaintype='G'
    smodel=[]
    calmode='p'
    append=False
    docallib=False
    gaintable=GainTables
    gainfield=['']
    interp=['']
    spwmap=[]
    parang=False
    async=False
    gaincal()


logprint ("Initial gain calibration on BP calibrator complete", logfileout='logs/semiFinalBPdcals1.log')
logprint ("Plotting initial phase gain calibration on BP calibrator", logfileout='logs/semiFinalBPdcals1.log')


for ii in range(nplots):
    filename='BPinitialgainphase'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='BPinitialgain.g'
    xaxis='time'
    yaxis='phase'
    poln=''
    field=''
    antenna=antPlot
    spw=''
    timerange=''
    subplot=311
    overplot=False
    clearpanel='Auto'
    iteration='antenna'
    plotrange=[0,0,-180,180]
    showflags=False
    plotsymbol='o-'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()


syscommand='rm -rf BPcal.b'
os.system(syscommand)

BPGainTables=copy.copy(priorcals)
BPGainTables.append('delay.k')
BPGainTables.append('BPinitialgain.g')

if (cal3C84_bp == True):
    default('bandpass')
    vis=ms_active
    caltable='BPcal.b'
    field=bandpass_field_select_string
    spw=''
    selectdata=True
    uvrange=uvrange3C84
    scan=bandpass_scan_select_string
    solint='inf'
    combine='scan'
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=5.0
    solnorm=False
    bandtype='B'
    fillgaps=0
    smodel=[]
    append=False
    docallib=False
    gaintable=BPGainTables
    gainfield=['']
    interp=['']
    spwmap=[]
    parang=False
    async=False
    bandpass()
else:
    default('bandpass')
    vis=ms_active
    caltable='BPcal.b'
    field=bandpass_field_select_string
    spw=''
    selectdata=True
    uvrange=''
    scan=bandpass_scan_select_string
    solint='inf'
    combine='scan'
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=5.0
    solnorm=False
    bandtype='B'
    fillgaps=0
    smodel=[]
    append=False
    docallib=False
    gaintable=BPGainTables
    gainfield=['']
    interp=['']
    spwmap=[]
    parang=False
    async=False
    bandpass()

logprint ("Bandpass calibration complete", logfileout='logs/semiFinalBPdcals1.log')
flaggedSolnResult=getCalFlaggedSoln('BPcal.b')

logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/semiFinalBPdcals1.log')
logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/semiFinalBPdcals1.log')


# Plot BP solutions

logprint ("Plotting bandpass solutions", logfileout='logs/semiFinalBPdcals1.log')

tb.open('BPcal.b')
dataVarCol = tb.getvarcol('CPARAM')
flagVarCol = tb.getvarcol('FLAG')
tb.close()
rowlist = dataVarCol.keys()
nrows = len(rowlist)
maxmaxamp = 0.0
maxmaxphase = 0.0
for rrow in rowlist:
    dataArr = dataVarCol[rrow]
    flagArr = flagVarCol[rrow]
    amps=np.abs(dataArr)
    phases=np.arctan2(np.imag(dataArr),np.real(dataArr))
    good=np.logical_not(flagArr)
    tmparr=amps[good]
    if (len(tmparr)>0):
        maxamp=np.max(amps[good])
        if (maxamp>maxmaxamp):
            maxmaxamp=maxamp
    tmparr=np.abs(phases[good])
    if (len(tmparr)>0):
        maxphase=np.max(np.abs(phases[good]))*180./pi
        if (maxphase>maxmaxphase):
            maxmaxphase=maxphase
ampplotmax=maxmaxamp
phaseplotmax=maxmaxphase

for ii in range(nplots):
    filename='BPcal_amp'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='BPcal.b'
    xaxis='freq'
    yaxis='amp'
    poln=''
    field=''
    antenna=antPlot
    spw=''
    timerange=''
    subplot=311
    overplot=False
    clearpanel='Auto'
    iteration='antenna'
    plotrange=[0,0,0,ampplotmax]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()

for ii in range(nplots):
    filename='BPcal_phase'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='BPcal.b'
    xaxis='freq'
    yaxis='phase'
    poln=''
    field=''
    antenna=antPlot
    spw=''
    timerange=''
    subplot=311
    overplot=False
    clearpanel='Auto'
    iteration='antenna'
    plotrange=[0,0,-phaseplotmax,phaseplotmax]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()


logprint ("Plotting complete", logfileout='logs/semiFinalBPdcals1.log')


logprint ("Applying semi-final delay and BP calibrations to all calibrators", logfileout='logs/semiFinalBPdcals1.log')

AllCalTables=copy.copy(priorcals)
AllCalTables.append('delay.k')
AllCalTables.append('BPcal.b')

ntables=len(AllCalTables)

default('applycal')
vis=ms_active
field=''
spw=''
selectdata=True
scan=calibrator_scan_select_string
gaintable=AllCalTables
interp=['']
spwmap=[]
parang=False
calwt=[False]*ntables
flagbackup=False
async=False
applycal()


logprint ("Plot calibrated calibrators to check for further flagging/RFI", logfileout='logs/semiFinalBPdcals1.log')

# NB: have to find a way to get plotms to reload data to show
# flagged result on second run

syscommand='rm -rf semifinalcalibratedcals1.png'
os.system(syscommand)

#clearstat()

default('plotms')
vis=ms_active
xaxis='freq'
yaxis='amp'
ydatacolumn='corrected'
selectdata=True
scan=calibrator_scan_select_string
correlation=corrstring
averagedata=True
avgtime='1e8s'
avgscan=False
transform=False
extendflag=False
iteraxis=''
coloraxis='antenna2'
plotrange=[]
title=''
xlabel=''
ylabel=''
showmajorgrid=False
showminorgrid=False
plotfile='semifinalcalibratedcals1.png'
interactive=False
overwrite=True
showgui=False
async=False
plotms()

# Calculate fractions of flagged solutions for final QA2

flaggedDelaySolns=getCalFlaggedSoln('delay.k')
flaggedBPSolns=getCalFlaggedSoln('BPcal.b')

if (flaggedDelaySolns['all']['total'] > 0):
    if (flaggedDelaySolns['antmedian']['fraction'] > critfrac):
        QA2_delay='Partial'
    else:
        QA2_delay='Pass'
else:
    QA2_delay='Fail'

logprint ("QA2_delay: "+QA2_delay, logfileout='logs/semiFinalBPdcals1.log')

if (flaggedBPSolns['all']['total'] > 0):
    if (flaggedBPSolns['antmedian']['fraction'] > 0.2):
        QA2_BP='Partial'
    else:
        QA2_BP='Pass'
else:
    QA2_BP='Fail'

logprint ("QA2_BP: "+QA2_BP, logfileout='logs/semiFinalBPdcals1.log')

if (QA2_delay=='Fail' or QA2_BP=='Fail'):
    QA2_semiFinalBPdcals1='Fail'
elif (QA2_delay=='Partial' or QA2_BP=='Partial'):
    QA2_semiFinalBPdcals1='Partial'

logprint ("QA2 score: "+QA2_semiFinalBPdcals1, logfileout='logs/semiFinalBPdcals1.log')
logprint ("Finished EVLA_pipe_semiFinalBPdcals1.py", logfileout='logs/semiFinalBPdcals1.log')
time_list=runtiming('semiFinalBPdcals1', 'end')

pipeline_save()
