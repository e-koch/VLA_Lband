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

logprint ("Starting EVLA_pipe_testBPdcals.py", logfileout='logs/testBPdcals.log')
time_list=runtiming('testBPdcals', 'start')
QA2_testBPdcals='Pass'

# INITIAL TEST CALIBRATIONS USING BANDPASS AND DELAY CALIBRATORS

logprint ("Finding a reference antenna", logfileout='logs/testBPdcals.log')

refantspw=''
refantfield=calibrator_field_select_string

findrefant=RefAntHeuristics(vis=ms_active,field=refantfield,geometry=True,flagging=True)
RefAntOutput=findrefant.calculate()
refAnt=str(RefAntOutput[0])

logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/testBPdcals.log')
logprint ("Doing test calibrations", logfileout='logs/testBPdcals.log')

# Do initial phase solutions on the delay calibrator

syscommand='rm -rf testdelayinitialgain.g'
os.system(syscommand)

if (cal3C84_d == True):
    default('gaincal')
    vis = ms_active
    caltable = 'testdelayinitialgain.g'
    field = delay_field_select_string
    spw = tst_delay_spw
    intent = ''
    selectdata = True
    uvrange = uvrange3C84
    scan = delay_scan_select_string
    solint = 'int'
    combine = 'scan'
    preavg = -1.0
    refant = refAnt
    minblperant = minBL_for_cal
    minsnr = 3.0
    solnorm = False
    gaintype = 'G'
    smodel = []
    calmode = 'p'
    append = False
    docallib = False
    gaintable = priorcals
    gainfield = ['']
    interp = ['']
    spwmap = []
    parang = False
    async = False
    gaincal()
else:
    default('gaincal')
    vis = ms_active
    caltable = 'testdelayinitialgain.g'
    field = delay_field_select_string
    spw = ''
    intent = ''
    selectdata = True
    uvrange = ''
    scan = delay_scan_select_string
    solint = 'int'
    combine = 'scan'
    preavg = -1.0
    refant = refAnt
    minblperant = minBL_for_cal
    minsnr = 3.0
    solnorm = False
    gaintype = 'G'
    smodel = []
    calmode = 'p'
    append = False
    docallib = False
    gaintable = priorcals
    gainfield = ['']
    interp = ['']
    spwmap = []
    parang = False
    async = False
    gaincal()


logprint ("Initial phase calibration on delay calibrator complete", logfileout='logs/testBPdcals.log')

# Do initial test delay calibration ("test" because more flagging may be
# needed for the final version)
# For the future: investigate multiband delay

syscommand='rm -rf testdelay.k'
os.system(syscommand)

flaggedSolnResult=testdelays(ms_active,'testdelay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/testBPdcals.log')
logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

if (flaggedSolnResult['all']['total'] > 0):
    fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
else:
    fracFlaggedSolns=1.0

# NB: in case the reference antenna has a bad baseband/IF, check
# a couple of reference antennas if there is a high fraction of
# flagged solutions

if (fracFlaggedSolns > critfrac):
    logprint ("Not enough good solutions, trying a different reference antenna", logfileout='logs/testBPdcals.log')
    refAnt=str(RefAntOutput[1])
    logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/testBPdcals.log')
    flaggedSolnResult=testdelays(ms_active,'testdelay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
    logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/testBPdcals.log')
    logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

    if (flaggedSolnResult['all']['total'] > 0):
        fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
    else:
        fracFlaggedSolns=1.0

    if (fracFlaggedSolns > critfrac):
        logprint ("Not enough good solutions, trying a different reference antenna", logfileout='logs/testBPdcals.log')
        refAnt=str(RefAntOutput[2])
        logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/testBPdcals.log')

        flaggedSolnResult=testdelays(ms_active,'testdelay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
        logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/testBPdcals.log')
        logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

        if (flaggedSolnResult['all']['total'] > 0):
            fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
        else:
            fracFlaggedSolns=1.0

        if (fracFlaggedSolns > critfrac):
            logprint ("Not enough good solutions, trying a different reference antenna", logfileout='logs/testBPdcals.log')
            refAnt=str(RefAntOutput[3])
            logprint ("The pipeline will use antenna "+refAnt+" as the reference", logfileout='logs/testBPdcals.log')
            flaggedSolnResult=testdelays(ms_active,'testdelay.k',delay_field_select_string,delay_scan_select_string,refAnt,minBL_for_cal,priorcals,cal3C84_d,uvrange3C84)
            logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/testBPdcals.log')
            logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

            if (flaggedSolnResult['all']['total'] > 0):
                fracFlaggedSolns=flaggedSolnResult['antmedian']['fraction']
            else:
                fracFlaggedSolns=1.0

            if (fracFlaggedSolns > critfrac):
                logprint ("WARNING, tried several reference antennas, there might be something wrong with your data", logfileout='logs/testBPdcals.log')

logprint ("Plotting test delays", logfileout='logs/testBPdcals.log')

nplots=int(numAntenna/3)

if ((numAntenna%3)>0):
    nplots = nplots + 1

for ii in range(nplots):
    filename='testdelay'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='testdelay.k'
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

# Do initial amplitude and phase gain solutions on the BPcalibrator and delay
# calibrator; the amplitudes are used for flagging; only phase
# calibration is applied in final BP calibration, so that solutions are
# not normalized per spw and take out the baseband filter shape

# Try running with solint of int_time, 3*int_time, and 10*int_time.
# If there is still a large fraction of failed solutions with
# solint=10*int_time the source may be too weak, and calibration via the
# pipeline has failed; will need to implement a mode to cope with weak
# calibrators (later)

if (delay_scan_select_string == bandpass_scan_select_string):
   testgainscans=bandpass_scan_select_string
else:
   testgainscans=bandpass_scan_select_string+','+delay_scan_select_string

if ((cal3C84_d == True) or (cal3C84_bp == True)):
    cal3C84=True
else:
    cal3C84=False

syscommand='rm -rf testBPdinitialgain.g'
os.system(syscommand)

soltime=int_time
solint='int'
flaggedSolnResult1=testBPdgains(ms_active,'testBPdinitialgain.g',tst_bpass_spw,testgainscans,solint,refAnt,minBL_for_cal,priorcals,cal3C84,uvrange3C84)
logprint("For solint = "+solint+" fraction of flagged solutions = "+str(flaggedSolnResult1['all']['fraction']), logfileout='logs/testBPdcals.log')
logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult1['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

if (flaggedSolnResult1['all']['total'] > 0):
    fracFlaggedSolns1=flaggedSolnResult1['antmedian']['fraction']
else:
    fracFlaggedSolns1=1.0

gain_solint1=solint
shortsol1=soltime

if (fracFlaggedSolns1 > 0.05):
    soltime=3.0*int_time
    solint=str(soltime)+'s'
    flaggedSolnResult3=testBPdgains(ms_active,'testBPdinitialgain3.g',tst_bpass_spw,testgainscans,solint,refAnt,minBL_for_cal,priorcals,cal3C84,uvrange3C84)
    logprint("For solint = "+solint+" fraction of flagged solutions = "+str(flaggedSolnResult3['all']['fraction']), logfileout='logs/testBPdcals.log')
    logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult3['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

    if (flaggedSolnResult3['all']['total'] > 0):
        fracFlaggedSolns3=flaggedSolnResult3['antmedian']['fraction']
    else:
        fracFlaggedSolns3=1.0

    if (fracFlaggedSolns3 < fracFlaggedSolns1):
        gain_solint1=solint
        shortsol1=soltime
        syscommand='rm -rf testBPdinitialgain.g'
        os.system(syscommand)
        syscommand='mv testBPdinitialgain3.g testBPdinitialgain.g'
        os.system(syscommand)
        if (fracFlaggedSolns3 > 0.05):
            soltime=10.0*int_time
            solint=str(soltime)+'s'
            flaggedSolnResult10=testBPdgains(ms_active,'testBPdinitialgain10.g',tst_bpass_spw,testgainscans,solint,refAnt,minBL_for_cal,priorcals,cal3C84,uvrange3C84)
            logprint("For solint = "+solint+" fraction of flagged solutions = "+str(flaggedSolnResult10['all']['fraction']), logfileout='logs/testBPdcals.log')
            logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult10['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

            if (flaggedSolnResult10['all']['total'] > 0):
                fracFlaggedSolns10=flaggedSolnResult10['antmedian']['fraction']
            else:
                fracFlaggedSolns10=1.0

            if (fracFlaggedSolns10 < fracFlaggedSolns3):
                gain_solint1=solint
                shortsol1=soltime
                syscommand='rm -rf testBPdinitialgain.g'
                os.system(syscommand)
                syscommand='mv testBPdinitialgain10.g testBPdinitialgain.g'
                os.system(syscommand)
                if (fracFlaggedSolns10 > 0.05):
                    logprint ("WARNING, large fraction of flagged solutions, there might be something wrong with your data", logfileout='logs/testBPdcals.log')

logprint ("Test amp and phase calibration on delay and bandpass calibrators complete", logfileout='logs/testBPdcals.log')
logprint ("Using short solint = "+gain_solint1, logfileout='logs/testBPdcals.log')

# Plot amplitude gain solutions

logprint ("Plotting amplitude gain solutions", logfileout='logs/testBPdcals.log')

tb.open('testBPdinitialgain.g')
cpar=tb.getcol('CPARAM')
flgs=tb.getcol('FLAG')
tb.close()
amps=np.abs(cpar)
good=np.logical_not(flgs)
maxamp=np.max(amps[good])
plotmax=maxamp

for ii in range(nplots):
    filename='testBPdinitialgainamp'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
#
    antPlot=str(ii*3)+'~'+str(ii*3+2)
#
    default('plotcal')
    caltable='testBPdinitialgain.g'
    xaxis='time'
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
    plotrange=[0,0,0,plotmax]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()


# Plot phase gain solutions

logprint ("Plotting phase gain solutions", logfileout='logs/testBPdcals.log')

for ii in range(nplots):
    filename='testBPdinitialgainphase'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
#
    antPlot=str(ii*3)+'~'+str(ii*3+2)
#
    default('plotcal')
    caltable='testBPdinitialgain.g'
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


# Now do test BPcal

logprint ("Doing test bandpass calibration", logfileout='logs/testBPdcals.log')

syscommand='rm -rf testBPcal.b'
os.system(syscommand)

BPGainTables=copy.copy(priorcals)
BPGainTables.append('testdelay.k')
BPGainTables.append('testBPdinitialgain.g')

if (cal3C84_bp == True):
    default('bandpass')
    vis=ms_active
    caltable='testBPcal.b'
    field=bandpass_field_select_string
    spw=''
    intent=''
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
    caltable='testBPcal.b'
    field=bandpass_field_select_string
    spw=''
    intent=''
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

logprint ("Test bandpass calibration complete", logfileout='logs/testBPdcals.log')
flaggedSolnResult=getCalFlaggedSoln('testBPcal.b')

logprint("Fraction of flagged solutions = "+str(flaggedSolnResult['all']['fraction']), logfileout='logs/testBPdcals.log')
logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

# Plot BP solutions and check for missing spws, antennas, etc.

tb.open('testBPcal.b')
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
    filename='testBPcal_amp'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
#
    antPlot=str(ii*3)+'~'+str(ii*3+2)
#
    default('plotcal')
    caltable='testBPcal.b'
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
    filename='testBPcal_phase'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
#
    antPlot=str(ii*3)+'~'+str(ii*3+2)
#
    default('plotcal')
    caltable='testBPcal.b'
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

logprint ("Plotting of test bandpass solutions complete", logfileout='logs/testBPdcals.log')

# Do blcal to take out closure errors from source structure for plotting
# NB: would be good to be able to specify smodel=[1,0,0,0] here since
# otherwise this only works for sources that are not primary calibrators
# NB: blcal crashes if a spw is missing from gaintable, can't use this
# for now

#default('blcal')
#vis=ms_active
#caltable='testBPblcal.bl'
#field=''
#spw=''
#selectdata=True
#scans=testgainscans
#solint='30min'
#combine='scan'
#freqdep=False
#calmode='ap'
#solnorm=False
#gaintable=[priorcals,'testdelay.k','testBPdinitialgain.g','testBPcal.b']
#gainfield=['']
#interp=['']
#spwmap=[]
#parang=False
#async=False
#blcal()

# NB: level of blcal corrections are an indicator of pointy-ness of
# BPcal and delay cal and/or system health

# Apply gain and bandpass solutions and inspect calibrated BP and delay
# calibrator data for RFI or other problems

logprint ("Applying test calibrations to BP and delay calibrators", logfileout='logs/testBPdcals.log')

AllCalTables=copy.copy(priorcals)
AllCalTables.append('testdelay.k')
AllCalTables.append('testBPdinitialgain.g')
AllCalTables.append('testBPcal.b')

ntables=len(AllCalTables)

default('applycal')
vis=ms_active
field=''
spw=''
intent=''
selectdata=True
scan=testgainscans
docallib=False
gaintable=AllCalTables
interp=['']
spwmap=[]
calwt=[False]*ntables
parang=False
flagbackup=False
async=False
applycal()

logprint ("Plot calibrated bandpass and delay calibrators", logfileout='logs/testBPdcals.log')

syscommand='rm -rf testcalibratedBPcal.png'
os.system(syscommand)

default('plotms')
vis=ms_active
xaxis='freq'
yaxis='amp'
ydatacolumn='corrected'
selectdata=True
field=bandpass_field_select_string
scan=bandpass_scan_select_string
correlation=corrstring
averagedata=True
avgtime='1e8s'
avgscan=True
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
plotfile='testcalibratedBPcal.png'
overwrite=True
showgui=False
async=False
plotms()

# Plot calibrated delay calibrator, if different from BP cal

if (delay_scan_select_string != bandpass_scan_select_string):
    syscommand='rm -rf testcalibrated_delaycal.png'
    os.system(syscommand)

    default('plotms')
    vis=ms_active
    xaxis='freq'
    yaxis='amp'
    ydatacolumn='corrected'
    selectdata=True
    scan=delay_scan_select_string
    correlation=corrstring
    averagedata=True
    avgtime='1e8s'
    avgscan=True
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
    plotfile='testcalibrated_delaycal.png'
    overwrite=True
    showgui=False
    async=False
    plotms()

# Calculate fractions of flagged solutions for final QA2

flaggedDelaySolns=getCalFlaggedSoln('testdelay.k')
flaggedGainSolns=getCalFlaggedSoln('testBPdinitialgain.g')
flaggedBPSolns=getCalFlaggedSoln('testBPcal.b')

if (flaggedDelaySolns['all']['total'] > 0):
    if (flaggedDelaySolns['antmedian']['fraction'] > critfrac):
        QA2_delay='Partial'
    else:
        QA2_delay='Pass'
else:
    QA2_delay='Fail'

logprint ("QA2_delay: "+QA2_delay, logfileout='logs/testBPdcals.log')

if (flaggedGainSolns['all']['total'] > 0):
    if (flaggedGainSolns['antmedian']['fraction'] > 0.1):
        QA2_gain='Partial'
    else:
        QA2_gain='Pass'
else:
    QA2_gain='Fail'

logprint ("QA2_gain: "+QA2_gain, logfileout='logs/testBPdcals.log')

if (flaggedBPSolns['all']['total'] > 0):
    if (flaggedBPSolns['antmedian']['fraction'] > 0.2):
        QA2_BP='Partial'
    else:
        QA2_BP='Pass'
else:
    QA2_BP='Fail'

logprint ("QA2_BP: "+QA2_BP, logfileout='logs/testBPdcals.log')

if (QA2_delay=='Fail' or QA2_gain=='Fail' or QA2_BP=='Fail'):
    QA2_testBPdcals='Fail'
elif (QA2_delay=='Partial' or QA2_gain=='Partial' or QA2_BP=='Partial'):
    QA2_testBPdcals='Partial'


logprint ("QA2 score: "+QA2_testBPdcals, logfileout='logs/testBPdcals.log')
logprint ("Finished EVLA_pipe_testBPdcals.py", logfileout='logs/testBPdcals.log')
time_list=runtiming('testBPdcals', 'end')

pipeline_save()
