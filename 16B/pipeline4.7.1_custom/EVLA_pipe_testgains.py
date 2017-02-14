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

# DO TEST GAIN CALIBRATIONS TO SEE IF MORE FLAGGING IS NEEDED AND TO
# ESTABLISH SHORT AND LONG SOLINTS
# (Needs some work to automate; note also that plotcal holds onto
# testgaincal.g in the table cache unless it has been exited using
# the gui, so only plot the final versions)

logprint ("Starting EVLA_pipe_testgains.py", logfileout='logs/testgains.log')
time_list=runtiming('testgains', 'start')
QA2_testgains='Pass'

print ""
print "Finding a reference antenna for gain calibrations"
print ""

refantspw=''
refantfield=calibrator_field_select_string

# NB: would use ms_active below instead of calibrators.ms when selection
# to exclude flagged data is implemented

findrefant=RefAntHeuristics(vis='calibrators.ms',field=refantfield,geometry=True,flagging=True)
RefAntOutput=findrefant.calculate()
refAnt=str(RefAntOutput[0])+','+str(RefAntOutput[1])+','+str(RefAntOutput[2])+','+str(RefAntOutput[3])

logprint ("The pipeline will use antenna(s) "+refAnt+" as the reference", logfileout='logs/testgains.log')

logprint ("Doing test gain calibration", logfileout='logs/testgains.log')

# First determine short solint for gain calibrator, and see if it is
# shorter or longer than gain_solint1 (determined on BPd cals)

# Start with solint='int'

syscommand='rm -rf testgaincal.g'
os.system(syscommand)

soltime=int_time
solint='int'
tst_gcal_spw=''
combtime='scan'

flaggedSolnResult1=testgains('calibrators.ms','testgaincal.g',tst_gcal_spw,calibrator_scan_select_string,solint,refAnt,minBL_for_cal,combtime)
logprint("For solint = "+solint+" fraction of flagged solutions = "+str(flaggedSolnResult1['all']['fraction']), logfileout='logs/testBPdcals.log')
logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult1['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

if (flaggedSolnResult1['all']['total'] > 0):
    fracFlaggedSolns1=flaggedSolnResult1['antmedian']['fraction']
else:
    fracFlaggedSolns1=1.0

shortsol2=soltime

if (fracFlaggedSolns1 > 0.05):
    soltime=3.0*int_time
    solint=str(soltime)+'s'
    flaggedSolnResult3=testgains('calibrators.ms','testgaincal3.g',tst_gcal_spw,calibrator_scan_select_string,solint,refAnt,minBL_for_cal,combtime)
    logprint("For solint = "+solint+" fraction of flagged solutions = "+str(flaggedSolnResult3['all']['fraction']), logfileout='logs/testBPdcals.log')
    logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult3['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

    if (flaggedSolnResult3['all']['total'] > 0):
        fracFlaggedSolns3=flaggedSolnResult3['antmedian']['fraction']
    else:
        fracFlaggedSolns3=1.0

    if (fracFlaggedSolns3 < fracFlaggedSolns1):
        shortsol2=soltime
        syscommand='rm -rf testgaincal.g'
        os.system(syscommand)
        syscommand='mv testgaincal3.g testgaincal.g'
        os.system(syscommand)
        if (fracFlaggedSolns3 > 0.05):
            soltime=10.0*int_time
            solint=str(soltime)+'s'
            flaggedSolnResult10=testgains('calibrators.ms','testgaincal10.g',tst_gcal_spw,calibrator_scan_select_string,solint,refAnt,minBL_for_cal,combtime)
            logprint("For solint = "+solint+" fraction of flagged solutions = "+str(flaggedSolnResult10['all']['fraction']), logfileout='logs/testBPdcals.log')
            logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResult10['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

            if (flaggedSolnResult10['all']['total'] > 0):
                fracFlaggedSolns10=flaggedSolnResult10['antmedian']['fraction']
            else:
                fracFlaggedSolns10=1.0

            if (fracFlaggedSolns10 < fracFlaggedSolns3):
                shortsol2=soltime
                syscommand='rm -rf testgaincal.g'
                os.system(syscommand)
                syscommand='mv testgaincal10.g testgaincal.g'
                os.system(syscommand)
                if (fracFlaggedSolns10 > 0.05):
                    solint='inf'
                    combtime=''
                    flaggedSolnResultScan=testgains('calibrators.ms','testgaincalscan.g',tst_gcal_spw,calibrator_scan_select_string,solint,refAnt,minBL_for_cal,combtime)
                    logprint("For solint = "+solint+" fraction of flagged solutions = "+str(flaggedSolnResultScan['all']['fraction']), logfileout='logs/testBPdcals.log')
                    logprint("Median fraction of flagged solutions per antenna = "+str(flaggedSolnResultScan['antmedian']['fraction']), logfileout='logs/testBPdcals.log')

                    if (flaggedSolnResultScan['all']['total'] > 0):
                        fracFlaggedSolnsScan=flaggedSolnResultScan['antmedian']['fraction']
                    else:
                        fracFlaggedSolnsScan=1.0

                    if (fracFlaggedSolnsScan < fracFlaggedSolns10):
                        shortsol2=longsolint
                        syscommand='rm -rf testgaincal.g'
                        os.system(syscommand)
                        syscommand='mv testgaincalscan.g testgaincal.g'
                        os.system(syscommand)
                        if (fracFlaggedSolnsScan > 0.05):
                            logprint ("Warning, large fraction of flagged solutions, there might be something wrong with your data", logfileout='logs/testBPdcals.log')

# determine max (shortsol1, shortsol2)

short_solint=max(shortsol1,shortsol2)
new_gain_solint1=str(short_solint)+'s'

logprint ("Using short solint = "+new_gain_solint1, logfileout='logs/testBPdcals.log')


# Plot solutions


logprint ("Plotting gain solutions", logfileout='logs/testgains.log')

nplots=int(numAntenna/3)

if ((numAntenna%3)>0):
    nplots = nplots + 1

tb.open('testgaincal.g')
cpar=tb.getcol('CPARAM')
flgs=tb.getcol('FLAG')
tb.close()
amps=np.abs(cpar)
good=np.logical_not(flgs)
maxamp=np.max(amps[good])
plotmax=maxamp

for ii in range(nplots):
    filename='testgaincal_amp'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
#
    antPlot=str(ii*3)+'~'+str(ii*3+2)
#
    default('plotcal')
    caltable='testgaincal.g'
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
    plotcal()

for ii in range(nplots):
    filename='testgaincal_phase'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)
#
    antPlot=str(ii*3)+'~'+str(ii*3+2)
#
    default('plotcal')
    caltable='testgaincal.g'
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
    plotcal()


logprint ("Plotting finished", logfileout='logs/testgains.log')


# Calculate fractions of flagged solutions for final QA2

flaggedGainSolns=getCalFlaggedSoln('testgaincal.g')

if (flaggedGainSolns['all']['total'] == 0):
    QA2_testgains='Fail'
elif (flaggedGainSolns['antmedian']['fraction'] > 0.1):
    QA2_testgains='Partial'

logprint ("QA2 score: "+QA2_testgains, logfileout='logs/testgains.log')
logprint ("Finished EVLA_pipe_testgains.py", logfileout='logs/testgains.log')
time_list=runtiming('testgains', 'end')

pipeline_save()

