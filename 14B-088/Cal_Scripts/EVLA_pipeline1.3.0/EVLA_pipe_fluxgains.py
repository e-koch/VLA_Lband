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

# MAKE GAIN TABLE FOR FLUX DENSITY BOOTSTRAPPING

# Make a gain table that includes gain and opacity corrections for final
# amp cal, for flux density bootstrapping

logprint ("Starting EVLA_pipe_fluxgains.py", logfileout='logs/fluxgains.log')
time_list=runtiming('fluxgains', 'start')
QA2_fluxgains='Pass'

#logprint ("Making fresh calibrators.ms", logfileout='logs/fluxgains.log')
#
#syscommand='rm -rf calibrators.ms'
#os.system(syscommand)
#
#default('split')
#vis=ms_active
#outputvis='calibrators.ms'
#datacolumn='corrected'
#field=''
#spw=''
#width=int(max(channels))
#antenna=''
#timebin='0s'
#timerange=''
#scan=calibrator_scan_select_string
#intent=''
#array=''
#uvrange=''
#correlation=''
#observation=''
#keepflags=False
#async=False
#split()

logprint ("Setting models for standard primary calibrators", logfileout='logs/fluxgains.log')

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
            logprint ("Center freq for spw "+str(myspw)+" = "+str(reference_frequency)+", observing band = "+EVLA_band, logfileout='logs/fluxgains.log')

            model_image = standard_source_names[ii]+'_'+EVLA_band+'.im'

            logprint ("Setting model for field "+str(myfield)+" spw "+str(myspw)+" using "+model_image, logfileout='logs/fluxgains.log')

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

logprint ("Making gain tables for flux density bootstrapping", logfileout='logs/fluxgains.log')

logprint ("Short solint = "+new_gain_solint1, logfileout='logs/fluxgains.log')
logprint ("Long solint = "+gain_solint2, logfileout='logs/fluxgains.log')

print ""
print "Finding a reference antenna"
print ""

refantspw=''
refantfield=calibrator_field_select_string

findrefant=RefAntHeuristics(vis='calibrators.ms',field=refantfield,geometry=True,flagging=True)
RefAntOutput=findrefant.calculate()
refAnt=str(RefAntOutput[0])+','+str(RefAntOutput[1])+','+str(RefAntOutput[2])+','+str(RefAntOutput[3])


logprint ("The pipeline will use antenna(s) "+refAnt+" as the reference", logfileout='logs/fluxgains.log')

# Derive amp gain table.  Note that gaincurves and opacity
# corrections have already been applied during applycal and split in
# semiFinalBPdcals/solint.py.

# Need to add check for 3C84 in here, when heuristics have been sorted out

default('gaincal')
vis='calibrators.ms'
caltable='fluxphaseshortgaincal.g'
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
#gaintable=filter(None, [priorcals,'delay.k','BPcal.b'])
gaintable=['']
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()

default('gaincal')
vis='calibrators.ms'
caltable='fluxgaincal.g'
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
#gaintable=filter(None, [priorcals,'delay.k','BPcal.b','fluxphaseshortgaincal.g'])
gaintable=['fluxphaseshortgaincal.g']
gainfield=['']
interp=['']
spwmap=[]
parang=False
async=False
gaincal()


logprint ("Gain table fluxgaincal.g is ready for flagging", logfileout='logs/fluxgains.log')

# Calculate fractions of flagged solutions for final QA2; note, can
# tolerate higher fraction of flagged solutions for this step than in
# other gain tables

flaggedGainSolns=getCalFlaggedSoln('fluxgaincal.g')

if (flaggedGainSolns['all']['total'] == 0):
    QA2_fluxgains='Fail'
elif (flaggedGainSolns['antmedian']['fraction'] > 0.2):
    QA2_fluxgains='Partial'

logprint ("QA2 score: "+QA2_fluxgains, logfileout='logs/fluxgains.log')
logprint ("Finished EVLA_pipe_fluxgains.py", logfileout='logs/fluxgains.log')
time_list=runtiming('fluxgains', 'end')

pipeline_save()
