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

# MAKE FINAL GAIN CALIBRATION TABLES

logprint ("Starting EVLA_pipe_finalcals.py", logfileout='logs/finalcals.log')
time_list=runtiming('finalcals', 'start')
QA2_finalcals='Pass'

# Find reference antenna again in case there has been more flagging

refantspw=''
refantfield=calibrator_field_select_string

findrefant=RefAntHeuristics(vis=ms_active,field=refantfield,geometry=True,flagging=True)
RefAntOutput=findrefant.calculate()
refAnt=str(RefAntOutput[0])+','+str(RefAntOutput[1])+','+str(RefAntOutput[2])+','+str(RefAntOutput[3])

logprint ("The pipeline will use antenna(s) "+refAnt+" as the reference", logfileout='logs/finalcals.log')


# Initial phase solutions on delay calibrator

syscommand='rm -rf finaldelayinitialgain.g'
os.system(syscommand)

if (cal3C84_d == True):
    default('gaincal')
    vis=ms_active
    caltable='finaldelayinitialgain.g'
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
    caltable='finaldelayinitialgain.g'
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

syscommand='rm -rf finaldelay.k'
os.system(syscommand)

GainTables=copy.copy(priorcals)
GainTables.append('finaldelayinitialgain.g')

if (cal3C84_d == True):
    default('gaincal')
    vis=ms_active
    caltable='finaldelay.k'
    field=''
    spw=''
    intent=''
    selectdata=True
    uvrange=uvrange3C84
    scan=delay_scan_select_string
    solint='inf'
    combine='scan'
    preavg=-1.0
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=3.0
    solnorm=False
    gaintype='K'
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
    caltable='finaldelay.k'
    field=''
    spw=''
    intent=''
    selectdata=True
    uvrange=''
    scan=delay_scan_select_string
    solint='inf'
    combine='scan'
    preavg=-1.0
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=3.0
    solnorm=False
    gaintype='K'
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


logprint ("Delay calibration complete", logfileout='logs/finalcals.log')


syscommand='rm -rf finalBPinitialgain.g'
os.system(syscommand)

GainTables=copy.copy(priorcals)
GainTables.append('finaldelay.k')

if (cal3C84_bp == True):
    default('gaincal')
    vis=ms_active
    caltable='finalBPinitialgain.g'
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
    caltable='finalBPinitialgain.g'
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



logprint ("Initial BP gain calibration complete", logfileout='logs/finalcals.log')


syscommand='rm -rf finalBPcal.b'
os.system(syscommand)

BPGainTables=copy.copy(priorcals)
BPGainTables.append('finaldelay.k')
BPGainTables.append('finalBPinitialgain.g')

if (cal3C84_bp == True):
    default('bandpass')
    vis=ms_active
    caltable='finalBPcal.b'
    field=bandpass_field_select_string
    spw=''
    selectdata=True
    uvrange=uvrange3C84
    scan=bandpass_scan_select_string
    solint='inf,4chan'
    combine='scan'
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=5.0
    solnorm=False
    bandtype='BPOLY'
    degamp=1
    degphase=1
    maskedge=20
    fillgaps=325
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
    caltable='finalBPcal.b'
    field=bandpass_field_select_string
    spw=''
    selectdata=True
    uvrange=''
    scan=bandpass_scan_select_string
    solint='inf,4chan'
    combine='scan'
    refant=refAnt
    minblperant=minBL_for_cal
    minsnr=5.0
    solnorm=False
    bandtype='BPOLY'
    degamp=1
    degphase=1
    maskedge=20
    fillgaps=325
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


logprint ("Bandpass calibration complete", logfileout='logs/finalcals.log')

AllCalTables=copy.copy(priorcals)
AllCalTables.append('finaldelay.k')
AllCalTables.append('finalBPcal.b')

#Derive an average phase solution for the bandpass calibrator to apply
#to all data to make QA plots easier to interpret.

default('gaincal')
vis=ms_active
caltable='averagephasegain.g'
field=bandpass_field_select_string
spw=''
selectdata=True
uvrange=''
scan=bandpass_scan_select_string
solint='inf'
combine='scan'
preavg=-1.0
refant=refAnt
minblperant=minBL_for_cal
minsnr=1.0
solnorm=False
gaintype='G'
smodel=[]
calmode='p'
append=False
docallib=False
gaintable=AllCalTables
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()

#In case any antenna is flagged by this process, unflag all solutions
#in this gain table (if an antenna does exist or has bad solutions from
#other steps, it will be flagged by those gain tables).

default('flagdata')
vis='averagephasegain.g'
mode='unflag'
action='apply'
flagbackup=False
savepars=False
async=False
flagdata()

AllCalTables=copy.copy(priorcals)
AllCalTables.append('finaldelay.k')
AllCalTables.append('finalBPcal.b')
AllCalTables.append('averagephasegain.g')

ntables=len(AllCalTables)

default('applycal')
vis=ms_active
field=''
spw=''
selectdata=True
scan=calibrator_scan_select_string
docallib=False
gaintable=AllCalTables
interp=['']
spwmap=[]
calwt=[False]*ntables
parang=False
flagbackup=False
async=False
applycal()


syscommand='rm -rf calibrators.ms'
os.system(syscommand)

default('split')
vis=ms_active
outputvis='calibrators.ms'
datacolumn='corrected'
field=''
spw=''
width=int(max(channels))
antenna=''
timebin='0s'
timerange=''
scan=calibrator_scan_select_string
intent=''
array=''
uvrange=''
correlation=''
observation=''
keepflags=False
async=False
split()

tb.open('calibrators.ms')

positions = []

for ii in range(0,len(field_positions[0][0])):
    positions.append([field_positions[0][0][ii], field_positions[1][0][ii]])

standard_source_names = [ '3C48', '3C138', '3C147', '3C286' ]
standard_source_fields = find_standards(positions)

ii=0
for fields in standard_source_fields:
    for myfield in fields:
        spws = field_spws[myfield]
        for myspw in spws:
            reference_frequency = center_frequencies[myspw]
            EVLA_band = find_EVLA_band(reference_frequency)
            logprint ("Center freq for spw "+str(myspw)+" = "+str(reference_frequency)+", observing band = "+EVLA_band, logfileout='logs/calprep.log')

            model_image = standard_source_names[ii]+'_'+EVLA_band+'.im'

            logprint ("Setting model for field "+str(myfield)+" spw "+str(myspw)+" using "+model_image, logfileout='logs/calprep.log')

            default('setjy')
            vis='calibrators.ms'
            field=str(myfield)
            spw=str(myspw)
            selectdata=False
            scalebychan=True
            standard='Perley-Butler 2013'
            model=model_image
            listmodels=False
            usescratch=False
            async=False
            setjy()
    ii=ii+1

tb.close()

fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

try:
    ff = open(fluxscale_output, 'r')
except IOError as err:
    logprint (fluxscale_output+" doesn't exist, error: "+err.filename, logfileout='logs/finalcals.log')


if calibrator_field_select_string == flux_field_select_string:
    logprint("No other calibrators found. No need to apply power-law fit to model column.",
             logfileout='logs/finalcals.log')

else:
    # looking for lines like:
    #2012-03-09 21:30:23     INFO    fluxscale::::    Flux density for J1717-3342 in SpW=3 is: 1.94158 +/- 0.0123058 (SNR = 157.777, N= 34)
    # sometimes they look like:
    #2012-03-09 21:30:23     INFO    fluxscale::::    Flux density for J1717-3342 in SpW=0 is:  INSUFFICIENT DATA
    # so watch for that.

    sources = []
    flux_densities = []
    spws = []

    #Find the field_ids in the dictionary returned from the CASA task fluxscale
    dictkeys = fluxscale_result.keys()
    keys_to_remove = ['freq', 'spwName', 'spwID']
    dictkeys = [field_id for field_id in dictkeys if field_id not in keys_to_remove]

    for field_id in dictkeys:
        sourcename = fluxscale_result[field_id]['fieldName']
        secondary_keys = fluxscale_result[field_id].keys()
        secondary_keys_to_remove=['fitRefFreq', 'spidxerr', 'spidx', 'fitFluxd', 'fieldName', 'fitFluxdErr']
        spwkeys = [spw_id for spw_id in secondary_keys if spw_id not in secondary_keys_to_remove]

        for spw_id in spwkeys:
            flux_d = list(fluxscale_result[field_id][spw_id]['fluxd'])
            flux_d_err = list(fluxscale_result[field_id][spw_id]['fluxdErr'])
            #spwslist  = list(int(spw_id))

            #flux_d = list(fluxscale_result[field_id]['fluxd'])
            #flux_d_err = list(fluxscale_result[field_id]['fluxdErr'])
            #spwslist  = list(fluxscale_result['spwID'])

            for i in range(0,len(flux_d)):
                if (flux_d[i] != -1.0 and flux_d[i] != 0.0):
                    sources.append(sourcename)
                    flux_densities.append([float(flux_d[i]), float(flux_d_err[i])])
                    spws.append(int(spw_id))

    ii = 0
    unique_sources = list(np.unique(sources))
    results = []
    for source in unique_sources:
        indices = []
        for ii in range(len(sources)):
            if (sources[ii] == source):
                indices.append(ii)
        bands = []
        for ii in range(len(indices)):
            bands.append(find_EVLA_band(center_frequencies[spws[indices[ii]]]))
        unique_bands = list(np.unique(bands))
        for band in unique_bands:
            lfreqs = []
            lfds = []
            lerrs = []
            uspws = []
            for ii in range(len(indices)):
                if find_EVLA_band(center_frequencies[spws[indices[ii]]]) == band:
                    lfreqs.append(log10(center_frequencies[spws[indices[ii]]]))
                    lfds.append(log10(flux_densities[indices[ii]][0]))
                    lerrs.append((flux_densities[indices[ii]][1])/(flux_densities[indices[ii]][0])/2.303)
                    uspws.append(spws[indices[ii]])

            if len(lfds) < 2:
               pfinal = [lfds[0], 0.0]
               covar = [0.0,0.0]
            else:
               alfds = scp.array(lfds)
               alerrs = scp.array(lerrs)
               alfreqs = scp.array(lfreqs)
               pinit = [0.0, 0.0]
               fit_out = scpo.leastsq(errfunc, pinit, args=(alfreqs, alfds, alerrs), full_output=1)
               pfinal = fit_out[0]
               covar = fit_out[1]
            aa = pfinal[0]
            bb = pfinal[1]
            reffreq = 10.0**lfreqs[0]/1.0e9
            fluxdensity = 10.0**(aa + bb*lfreqs[0])
            spix = bb
            results.append([ source, uspws, fluxdensity, spix, reffreq ])
            logprint(source + ' ' + band + ' fitted spectral index = ' + str(spix), logfileout='logs/finalcals.log')
            logprint("Frequency, data, and fitted data:", logfileout='logs/finalcals.log')
            for ii in range(len(lfreqs)):
                SS = fluxdensity * (10.0**lfreqs[ii]/reffreq/1.0e9)**spix
                logprint('    '+str(10.0**lfreqs[ii]/1.0e9)+'  '+ str(10.0**lfds[ii])+'  '+str(SS), logfileout='logs/finalcals.log')


    logprint ("Setting power-law fit in the model column", logfileout='logs/finalcals.log')


    for result in results:
        for spw_i in result[1]:
            logprint('Running setjy on spw '+str(spw_i), logfileout='logs/finalcals.log')
            default('setjy')
            vis='calibrators.ms'
            field = str(result[0])
            #spw = ','.join(["%s" % ii for ii in result[1]])
            spw = str(spw_i)
            selectdata=False
            scalebychan=True
            standard='manual'
            fluxdensity = [ result[2], 0, 0, 0 ]
            spix = result[3]
            reffreq = str(result[4])+'GHz'
            usescratch=False
            async=False
            setjy()

# Derive gain tables.  Note that gaincurves, opacity corrections and
# antenna position corrections have already been applied during applycal
# and split in above.

# Need to add check for 3C84 in here, when heuristics have been sorted out

default('gaincal')
vis='calibrators.ms'
caltable='phaseshortgaincal.g'
field=''
spw=''
intent=''
selectdata=False
solint=new_gain_solint1
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
#gaintable=filter(None, [priorcals,'finaldelay.k','finalBPcal.b'])
gaintable=['']
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()

default('gaincal')
vis='calibrators.ms'
caltable='finalampgaincal.g'
field=''
spw=''
intent=''
selectdata=False
solint=gain_solint2
combine='scan'
preavg=-1.0
refant=refAnt
minblperant=minBL_for_cal
minsnr=5.0
solnorm=False
gaintype='G'
smodel=[]
calmode='ap'
append=False
docallib=False
#gaintable=filter(None, [priorcals,'finaldelay.k','finalBPcal.b','phaseshortgaincal.g'])
gaintable=['phaseshortgaincal.g']
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()

default('gaincal')
vis='calibrators.ms'
caltable='finalphasegaincal.g'
field=''
spw=''
intent=''
selectdata=False
solint=gain_solint2
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
#gaintable=filter(None, [priorcals,'finaldelay.k','finalBPcal.b','finalampgaincal.g'])
gaintable=['finalampgaincal.g']
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()


logprint ("Final calibration tables created", logfileout='logs/finalcals.log')
logprint ("Plotting final calibration tables", logfileout='logs/finalcals.log')

# do some plotting

nplots=int(numAntenna/3)

if ((numAntenna%3)>0):
    nplots = nplots + 1

tb.open('finaldelay.k')
fpar = tb.getcol('FPARAM')
tb.close()
delays = np.abs(fpar)
maxdelay = np.max(delays)

for ii in range(nplots):
    filename='finaldelay'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finaldelay.k'
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

for ii in range(nplots):
    filename='finalBPinitialgainphase'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finalBPinitialgain.g'
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

# import numpy as np

# def bpoly_model(freqs, params):

#     split = len(params)/2

#     bpoly_curve = np.zeros((2, len(freqs)))
#     for i in range(split):
#         bpoly_curve[0, :] += params[i] * np.power(freqs, i)
#         bpoly_curve[1, :] += params[split+i] * np.power(freqs, i)

#     return bpoly_curve

# tb.open(ms_active+"/SPECTRAL_WINDOW")
# freqs = tb.getvarcol("CHAN_FREQ")
# tb.close()

# tb.open('finalBPcal.b')
# phaseVarCol = tb.getvarcol('POLY_COEFF_PHASE')
# ampVarCol = tb.getvarcol('POLY_COEFF_AMP')
# flagVarCol = tb.getvarcol('FLAG')
# tb.close()

# rowlist = ampVarCol.keys()
# maxmaxamp = 0.0
# maxmaxphase = 0.0
# nspw = 0
# for i, rrow in enumerate(rowlist):
#     if i == (nspw+1)*numAntenna:
#         nspw += 1
#     # Check if it's flagged
#     if not flagVarCol[rrow]:
#         continue
#     maxamp = np.max(bpoly_model(freqs[nspw], ampVarCol[rrow]))
#     maxphase = np.max(bpoly_model(freqs[nspw], phaseVarCol[rrow])) * 180./np.pi

#     if maxamp > maxmaxamp:
#         maxmaxamp = maxamp
#     if maxphase > maxmaxphase:
#         maxmaxphase = maxphase
# ampplotmax = maxmaxamp
# phaseplotmax = maxmaxphase

for ii in range(nplots):
    filename='finalBPcal_amp'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finalBPcal.b'
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
    # plotrange=[0,0,0,ampplotmax]
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
    filename='finalBPcal_phase'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finalBPcal.b'
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
    # plotrange=[0,0,-phaseplotmax,phaseplotmax]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()

# Plot the bandpasses per SPW as well

tb.open(ms_active+"/SPECTRAL_WINDOW")
nspws = tb.getcol("NAME").shape[0]
tb.close()

for ii in range(nspws):
    filename='finalBPcal_amp_spw_'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    default('plotcal')
    caltable='finalBPcal.b'
    xaxis='freq'
    yaxis='amp'
    poln=''
    field=''
    antenna=''
    spw=str(ii)
    timerange=''
    subplot=111
    overplot=False
    clearpanel='Auto'
    iteration=''
    # plotrange=[0,0,0,ampplotmax]
    showflags=False
    plotsymbol='o'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()

for ii in range(nspws):
    filename='finalBPcal_phase_spw_'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finalBPcal.b'
    xaxis='freq'
    yaxis='phase'
    poln=''
    field=''
    antenna=''
    spw=str(ii)
    timerange=''
    subplot=111
    overplot=False
    clearpanel='Auto'
    iteration=''
    # plotrange=[0,0,-phaseplotmax,phaseplotmax]
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
    filename='phaseshortgaincal'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='phaseshortgaincal.g'
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

tb.open('finalampgaincal.g')
cpar=tb.getcol('CPARAM')
flgs=tb.getcol('FLAG')
tb.close()
amps=np.abs(cpar)
good=np.logical_not(flgs)
maxamp=np.max(amps[good])
plotmax=max(2.0,maxamp)

for ii in range(nplots):
    filename='finalamptimecal'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finalampgaincal.g'
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
    plotsymbol='o-'
    plotcolor='blue'
    markersize=5.0
    fontsize=10.0
    showgui=False
    figfile=filename
    async=False
    plotcal()

for ii in range(nplots):
    filename='finalampfreqcal'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finalampgaincal.g'
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

for ii in range(nplots):
    filename='finalphasegaincal'+str(ii)+'.png'
    syscommand='rm -rf '+filename
    os.system(syscommand)

    antPlot=str(ii*3)+'~'+str(ii*3+2)

    default('plotcal')
    caltable='finalphasegaincal.g'
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


# Calculate fractions of flagged solutions for final QA2

flaggedDelaySolns=getCalFlaggedSoln('finaldelay.k')
# flaggedBPSolns=getCalFlaggedSoln('finalBPcal.b')
flaggedAmpSolns=getCalFlaggedSoln('finalampgaincal.g')
flaggedPhaseSolns=getCalFlaggedSoln('finalphasegaincal.g')

if (flaggedDelaySolns['all']['total'] > 0):
    if (flaggedDelaySolns['antmedian']['fraction'] > critfrac):
        QA2_delay='Partial'
    else:
        QA2_delay='Pass'
else:
    QA2_delay='Fail'

# if (flaggedBPSolns['all']['total'] > 0):
#     if (flaggedBPSolns['antmedian']['fraction'] > 0.2):
#         QA2_BP='Partial'
#     else:
#         QA2_BP='Pass'
# else:
#     QA2_BP='Fail'
QA2_BP = 'Pass'

if (flaggedAmpSolns['all']['total'] > 0):
    if (flaggedAmpSolns['antmedian']['fraction'] > 0.1):
        QA2_amp='Partial'
    else:
        QA2_amp='Pass'
else:
    QA2_amp='Fail'

if (flaggedPhaseSolns['all']['total'] > 0):
    if (flaggedPhaseSolns['antmedian']['fraction'] > 0.1):
        QA2_phase='Partial'
    else:
        QA2_phase='Pass'
else:
    QA2_phase='Fail'

if (QA2_delay=='Fail' or QA2_BP=='Fail' or QA2_amp=='Fail' or QA2_phase=='Fail'):
    QA2_finalcals='Fail'
elif (QA2_delay=='Partial' or QA2_BP=='Partial' or QA2_amp=='Partial' or QA2_phase=='Partial'):
    QA2_finalcals='Partial'

logprint ("QA2 score: "+QA2_finalcals, logfileout='logs/finalcals.log')
logprint ("Finished EVLA_pipe_finalcals.py", logfileout='logs/finalcals.log')
time_list=runtiming('finalcals', 'end')

pipeline_save()
