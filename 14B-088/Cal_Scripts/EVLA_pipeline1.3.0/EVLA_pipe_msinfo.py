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
import casa

logprint ("Starting EVLA_pipe_msinfo.py", logfileout='logs/msinfo.log')
time_list=runtiming('msinfo', 'start')
QA2_msinfo='Pass'

# Run listobs

logprint ("Listing ms contents", logfileout='logs/msinfo.log')

listname=msname.rstrip('ms') + 'listobs'
syscommand='rm -rf '+listname
os.system(syscommand)

default('listobs')
vis=ms_active
selectdata=False
verbose=True
listfile=listname
async=False
listobs()

# GET SOME INFORMATION FROM THE MS THAT WILL BE NEEDED LATER

# Identify spw information

tb.open(ms_active+'/SPECTRAL_WINDOW')
channels = tb.getcol('NUM_CHAN')
originalBBClist = tb.getcol('BBC_NO')
spw_bandwidths=tb.getcol('TOTAL_BANDWIDTH')
reference_frequencies = tb.getcol('REF_FREQUENCY')
center_frequencies = []
for ii in range(len(reference_frequencies)):
    center_frequencies.append(reference_frequencies[ii]+spw_bandwidths[ii]/2)
tb.close()
bands = []
for center_frequency in center_frequencies:
    bands.append(find_EVLA_band(center_frequency))
unique_bands = uniq(bands)
unique_bands_string = ','.join(["%s" % ii for ii in unique_bands])
logprint("unique band string = " + unique_bands_string, logfileout='logs/msinfo.log')

numSpws = len(channels)

# Set up spw selection for initial gain solutions

tst_delay_spw=''
all_spw=''

for ispw in range(numSpws):
    endch1=int(channels[ispw]/3.0)
    endch2=int(2.0*channels[ispw]/3.0)+1
    if (ispw<max(range(numSpws))):
        tst_delay_spw=tst_delay_spw+str(ispw)+':'+str(endch1)+'~'+str(endch2)+','
        all_spw=all_spw+str(ispw)+','
    else:
        tst_delay_spw=tst_delay_spw+str(ispw)+':'+str(endch1)+'~'+str(endch2)
        all_spw=all_spw+str(ispw)

tst_bpass_spw=tst_delay_spw

# Identify number of fields, positions, and source IDs

tb.open(ms_active+'/FIELD')
numFields = tb.nrows()
field_positions = tb.getcol('PHASE_DIR')
field_ids=range(numFields)
field_names=tb.getcol('NAME')
tb.close()

# Map field IDs to spws

field_spws = []
for ii in range(numFields):
    field_spws.append(spwsforfield(ms_active,ii))

# Identify scan numbers, map scans to field ID, and run scan summary
# (needed for figuring out integration time later)

tb.open(ms_active)
scanNums = sorted(np.unique(tb.getcol('SCAN_NUMBER')))
field_scans = []
for ii in range(0,numFields):
    subtable = tb.query('FIELD_ID==%s'%ii)
    field_scans.append(list(np.unique(subtable.getcol('SCAN_NUMBER'))))
tb.close()

## field_scans is now a list of lists containing the scans for each field.
## so, to access all the scans for the fields, you'd:
#
#for ii in range(0,len(field_scans)):
#   for jj in range(0,len(field_scans[ii]))
#
## the jj'th scan of the ii'th field is in field_scans[ii][jj]

# Identify intents

tb.open(ms_active+'/STATE')
intents=tb.getcol('OBS_MODE')
tb.close()

# Figure out integration time used

ms.open(ms_active)
scan_summary = ms.getscansummary()
ms_summary = ms.summary()
ms.close()
startdate=float(ms_summary['BeginTime'])
#
# scan list
#
integ_scan_list = []
for scan in scan_summary:
    integ_scan_list.append(int(scan))
sorted_scan_list = sorted(integ_scan_list)
#
# find max and median integration times
#
integration_times = []
for ii in sorted_scan_list:
    integration_times.append(scan_summary[str(ii)]['0']['IntegrationTime'])

maximum_integration_time = max(integration_times)
median_integration_time = np.median(integration_times)

if (maximum_integration_time != median_integration_time):
    logprint ("Warning:", logfileout='logs/msinfo.log')
    logprint ("Median integration time = "+str(median_integration_time), logfileout='logs/msinfo.log')
    logprint ("Maximum integration time = "+str(maximum_integration_time), logfileout='logs/msinfo.log')

int_time=maximum_integration_time
logprint ("Maximum integration time is "+str(int_time)+"s", logfileout='logs/msinfo.log')

# Find scans for quacking

scan_list = [1]
old_scan = scan_summary[str(sorted_scan_list[0])]['0']
old_field = old_scan['FieldId']
old_spws = old_scan['SpwIds']
for ii in range(1,len(sorted_scan_list)):
    new_scan = scan_summary[str(sorted_scan_list[ii])]['0']
    new_field = new_scan['FieldId']
    new_spws = new_scan['SpwIds']
    if ((new_field != old_field) or (set(new_spws) != set(old_spws))):
        scan_list.append(sorted_scan_list[ii])
        old_field = new_field
        old_spws = new_spws
quack_scan_string = ','.join(["%s" % ii for ii in scan_list])

# For 1 GHz wide basebands, figure out which spws are associated
# with the edges of the baseband filters

sorted_frequencies = sorted(reference_frequencies)
sorted_indices = []

for ii in range (0,len(sorted_frequencies)):
    for jj in range (0,len(reference_frequencies)):
        if (sorted_frequencies[ii] == reference_frequencies[jj]):
            sorted_indices.append(jj)

spwList = []
BBC_bandwidths = []
ii = 0

while (ii < len(sorted_frequencies)):
    upper_frequency = sorted_frequencies[ii] + spw_bandwidths[sorted_indices[ii]]
    BBC_bandwidth = spw_bandwidths[sorted_indices[ii]]
    thisSpwList = [sorted_indices[ii]]
    jj = ii + 1
    while (jj < len(sorted_frequencies)):
        lower_frequency = sorted_frequencies[jj]
        if ((fabs(lower_frequency - upper_frequency) < 1.0) and \
            (originalBBClist[sorted_indices[ii]] == originalBBClist[sorted_indices[jj]])):
            thisSpwList.append(sorted_indices[jj])
            upper_frequency += spw_bandwidths[sorted_indices[jj]]
            BBC_bandwidth += spw_bandwidths[sorted_indices[jj]]
            jj += 1
            ii += 1
        else:
            jj = len(sorted_frequencies)
    spwList.append(thisSpwList)
    BBC_bandwidths.append(BBC_bandwidth)
    ii += 1

# spwList is now a list of lists of contiguous spws, which have
# bandwidths in BBC_bandwidths

low_spws = []
high_spws = []

for ii in range(0,len(BBC_bandwidths)):
    if (BBC_bandwidths[ii] > 1.0e9):
        low_spws.append(spwList[ii][0])
        high_spws.append(spwList[ii][len(spwList[ii])-1])

logprint ("Bottom ends of baseband filters are spws: "+str(low_spws), logfileout='logs/msinfo.log')
logprint ("Top ends of baseband filters are spws: "+str(high_spws), logfileout='logs/msinfo.log')

if (len(low_spws) != len(high_spws)):
    logprint("Error! Something is wrong with the spw identification",logfileout='logs/msinfo.log')
    raise Exception("Error! Something is wrong with the spw identification")

# Identify scans and fields associated with different calibrator intents

# NB: the scan intent definitions changed in the OPT on Feb 21,
# 2012.  So test on date:

if startdate <= 55978.50:
    bandpass_state_IDs = []
    delay_state_IDs = []
    flux_state_IDs = []
    polarization_state_IDs = []
    phase_state_IDs = []
    calibrator_state_IDs = []
    pointing_state_IDs = []
    for state_ID in range(0,len(intents)):
        state_intents = intents[state_ID].rsplit(',')
        for intent in range(0,len(state_intents)):
            scan_intent = state_intents[intent].rsplit('#')[0]
            subscan_intent = state_intents[intent].rsplit('#')[1]
            if (scan_intent == 'CALIBRATE_BANDPASS'):
                bandpass_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_DELAY'):
                delay_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_AMPLI'):
                flux_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_POLARIZATION'):
                polarization_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_PHASE'):
                phase_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_POINTING'):
                pointing_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)

    tb.open(ms_active)

    if (len(flux_state_IDs) == 0):
        QA2_msinfo='Fail'
        logprint("ERROR: No flux density calibration scans found", logfileout='logs/msinfo.log')
        raise Exception("No flux density calibration scans found")
    else:
        flux_state_select_string = ('STATE_ID in [%s'%flux_state_IDs[0])
        for state_ID in range(1,len(flux_state_IDs)):
            flux_state_select_string += (',%s')%flux_state_IDs[state_ID]
        flux_state_select_string += ']'
        subtable = tb.query(flux_state_select_string)
        flux_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        flux_scan_select_string = ','.join(["%s" % ii for ii in flux_scan_list])
        logprint ("Flux density calibrator(s) scans are "+flux_scan_select_string, logfileout='logs/msinfo.log')
        flux_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        flux_field_select_string = ','.join(["%s" % ii for ii in flux_field_list])
        logprint ("Flux density calibrator(s) are fields "+flux_field_select_string, logfileout='logs/msinfo.log')

    if (len(bandpass_state_IDs) == 0):
        logprint ("No bandpass calibration scans defined, using flux density calibrator(s)")
        bandpass_scan_select_string=flux_scan_select_string
        logprint ("Bandpass calibrator(s) scans are "+bandpass_scan_select_string, logfileout='logs/msinfo.log')
        bandpass_field_select_string=flux_field_select_string
        logprint ("Bandpass calibrator(s) are fields "+bandpass_field_select_string, logfileout='logs/msinfo.log')
    else:
        bandpass_state_select_string = ('STATE_ID in [%s'%bandpass_state_IDs[0])
        for state_ID in range(1,len(bandpass_state_IDs)):
            bandpass_state_select_string += (',%s')%bandpass_state_IDs[state_ID]
        bandpass_state_select_string += ']'
        subtable = tb.query(bandpass_state_select_string)
        bandpass_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        bandpass_scan_select_string = ','.join(["%s" % ii for ii in bandpass_scan_list])
        logprint ("Bandpass calibrator(s) scans are "+bandpass_scan_select_string, logfileout='logs/msinfo.log')
        bandpass_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        bandpass_field_select_string = ','.join(["%s" % ii for ii in bandpass_field_list])
        logprint ("Bandpass calibrator(s) are fields "+bandpass_field_select_string, logfileout='logs/msinfo.log')
        if (len(bandpass_field_list) > 1):
            logprint ("WARNING: More than one field is defined as the bandpass calibrator.", logfileout='logs/msinfo.log')
            logprint ("WARNING: Models are required for all BP calibrators if multiple fields", logfileout='logs/msinfo.log')
            logprint ("WARNING: are to be used, not yet implemented; the pipeline will use", logfileout='logs/msinfo.log')
            logprint ("WARNING: only the first field.", logfileout='logs/msinfo.log')
            bandpass_field_select_string = str(bandpass_field_list[0])

    if (len(delay_state_IDs) == 0):
        logprint ("No delay calibration scans defined, using bandpass calibrator")
        delay_scan_select_string=bandpass_scan_select_string
        logprint ("Delay calibrator(s) scans are "+delay_scan_select_string, logfileout='logs/msinfo.log')
        delay_field_select_string=bandpass_field_select_string
        logprint ("Delay calibrator(s) are fields "+delay_field_select_string, logfileout='logs/msinfo.log')
    else:
        delay_state_select_string = ('STATE_ID in [%s'%delay_state_IDs[0])
        for state_ID in range(1,len(delay_state_IDs)):
            delay_state_select_string += (',%s')%delay_state_IDs[state_ID]
        delay_state_select_string += ']'
        subtable = tb.query(delay_state_select_string)
        delay_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        delay_scan_select_string = ','.join(["%s" % ii for ii in delay_scan_list])
        logprint ("Delay calibrator(s) scans are "+delay_scan_select_string, logfileout='logs/msinfo.log')
        delay_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        delay_field_select_string = ','.join(["%s" % ii for ii in delay_field_list])
        logprint ("Delay calibrator(s) are fields "+delay_field_select_string, logfileout='logs/msinfo.log')
        if (len(delay_field_list) > 1):
            logprint ("WARNING: More than one field is defined as the delay calibrator.", logfileout='logs/msinfo.log')
            logprint ("WARNING: Models are required for all delay calibrators if multiple fields", logfileout='logs/msinfo.log')
            logprint ("WARNING: are to be used, not yet implemented; the pipeline will use", logfileout='logs/msinfo.log')
            logprint ("WARNING: only the first field.", logfileout='logs/msinfo.log')
            delay_field_select_string = str(delay_field_list[0])

    if (len(polarization_state_IDs) == 0):
        logprint("No polarization calibration scans defined, no polarization calibration possible", logfileout='logs/msinfo.log')
        polarization_scan_select_string=''
        polarization_field_select_string=''
    else:
        logprint("Warning: polarization calibration scans found, but polarization calibration not yet implemented", logfileout='logs/msinfo.log')
        polarization_state_select_string = ('STATE_ID in [%s'%polarization_state_IDs[0])
        for state_ID in range(1,len(polarization_state_IDs)):
            polarization_state_select_string += (',%s')%polarization_state_IDs[state_ID]
        polarization_state_select_string += ']'
        subtable = tb.query(polarization_state_select_string)
        polarization_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        polarization_scan_select_string = ','.join(["%s" % ii for ii in polarization_scan_list])
        logprint ("Polarization calibrator(s) scans are "+polarization_scan_select_string, logfileout='logs/msinfo.log')
        polarization_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        polarization_field_select_string = ','.join(["%s" % ii for ii in polarization_field_list])
        logprint ("Polarization calibrator(s) are fields "+polarization_field_select_string, logfileout='logs/msinfo.log')

    if (len(phase_state_IDs) == 0):
        QA2_msinfo='Fail'
        logprint("ERROR: No gain calibration scans found", logfileout='logs/msinfo.log')
        raise Exception("No gain calibration scans found")
    else:
        phase_state_select_string = ('STATE_ID in [%s'%phase_state_IDs[0])
        for state_ID in range(1,len(phase_state_IDs)):
            phase_state_select_string += (',%s')%phase_state_IDs[state_ID]
        phase_state_select_string += ']'
        subtable = tb.query(phase_state_select_string)
        phase_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        phase_scan_select_string = ','.join(["%s" % ii for ii in phase_scan_list])
        logprint ("Phase calibrator(s) scans are "+phase_scan_select_string, logfileout='logs/msinfo.log')
        phase_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        phase_field_select_string = ','.join(["%s" % ii for ii in phase_field_list])
        logprint ("Phase calibrator(s) are fields "+phase_field_select_string, logfileout='logs/msinfo.log')

# Find all calibrator scans and fields

    calibrator_state_select_string = ('STATE_ID in [%s'%calibrator_state_IDs[0])
    for state_ID in range(1,len(calibrator_state_IDs)):
        calibrator_state_select_string += (',%s')%calibrator_state_IDs[state_ID]

    calibrator_state_select_string += ']'
    subtable = tb.query(calibrator_state_select_string)
    calibrator_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
    calibrator_scan_select_string = ','.join(["%s" % ii for ii in calibrator_scan_list])
    calibrator_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
    calibrator_field_select_string = ','.join(["%s" % ii for ii in calibrator_field_list])

    tb.close()
else:
    bandpass_state_IDs = []
    delay_state_IDs = []
    flux_state_IDs = []
    polarization_state_IDs = []
    phase_state_IDs = []
    amp_state_IDs = []
    calibrator_state_IDs = []
    pointing_state_IDs = []
    for state_ID in range(0,len(intents)):
        state_intents = intents[state_ID].rsplit(',')
        for intent in range(0,len(state_intents)):
            scan_intent = state_intents[intent].rsplit('#')[0]
            subscan_intent = state_intents[intent].rsplit('#')[1]
            if (scan_intent == 'CALIBRATE_BANDPASS'):
                bandpass_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_DELAY'):
                delay_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_FLUX'):
                flux_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_POLARIZATION'):
                polarization_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_AMPLI'):
                amp_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_PHASE'):
                phase_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)
            elif (scan_intent == 'CALIBRATE_POINTING'):
                pointing_state_IDs.append(state_ID)
                calibrator_state_IDs.append(state_ID)

    tb.open(ms_active)

    if (len(flux_state_IDs) == 0):
        QA2_msinfo='Fail'
        logprint("ERROR: No flux density calibration scans found", logfileout='logs/msinfo.log')
        raise Exception("No flux density calibration scans found")
    else:
        flux_state_select_string = ('STATE_ID in [%s'%flux_state_IDs[0])
        for state_ID in range(1,len(flux_state_IDs)):
            flux_state_select_string += (',%s')%flux_state_IDs[state_ID]
        flux_state_select_string += ']'
        subtable = tb.query(flux_state_select_string)
        flux_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        flux_scan_select_string = ','.join(["%s" % ii for ii in flux_scan_list])
        logprint ("Flux density calibrator(s) scans are "+flux_scan_select_string, logfileout='logs/msinfo.log')
        flux_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        flux_field_select_string = ','.join(["%s" % ii for ii in flux_field_list])
        logprint ("Flux density calibrator(s) are fields "+flux_field_select_string, logfileout='logs/msinfo.log')

    if (len(bandpass_state_IDs) == 0):
        logprint ("No bandpass calibration scans defined, using flux density calibrator", logfileout='logs/msinfo.log')
        bandpass_scan_select_string=flux_scan_select_string
        logprint ("Bandpass calibrator(s) scans are "+bandpass_scan_select_string, logfileout='logs/msinfo.log')
        bandpass_field_select_string=flux_field_select_string
        logprint ("Bandpass calibrator(s) are fields "+bandpass_field_select_string, logfileout='logs/msinfo.log')
    else:
        bandpass_state_select_string = ('STATE_ID in [%s'%bandpass_state_IDs[0])
        for state_ID in range(1,len(bandpass_state_IDs)):
            bandpass_state_select_string += (',%s')%bandpass_state_IDs[state_ID]
        bandpass_state_select_string += ']'
        subtable = tb.query(bandpass_state_select_string)
        bandpass_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        bandpass_scan_select_string = ','.join(["%s" % ii for ii in bandpass_scan_list])
        logprint ("Bandpass calibrator(s) scans are "+bandpass_scan_select_string, logfileout='logs/msinfo.log')
        bandpass_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        bandpass_field_select_string = ','.join(["%s" % ii for ii in bandpass_field_list])
        logprint ("Bandpass calibrator(s) are fields "+bandpass_field_select_string, logfileout='logs/msinfo.log')
        if (len(bandpass_field_list) > 1):
            logprint ("WARNING: More than one field is defined as the bandpass calibrator.", logfileout='logs/msinfo.log')
            logprint ("WARNING: Models are required for all BP calibrators if multiple fields", logfileout='logs/msinfo.log')
            logprint ("WARNING: are to be used, not yet implemented; the pipeline will use", logfileout='logs/msinfo.log')
            logprint ("WARNING: only the first field.", logfileout='logs/msinfo.log')
            bandpass_field_select_string = str(bandpass_field_list[0])

    if (len(delay_state_IDs) == 0):
        logprint ("No delay calibration scans defined, using bandpass calibrator", logfileout='logs/msinfo.log')
        delay_scan_select_string=bandpass_scan_select_string
        logprint ("Delay calibrator(s) scans are "+delay_scan_select_string, logfileout='logs/msinfo.log')
        delay_field_select_string=bandpass_field_select_string
        logprint ("Delay calibrator(s) are fields "+delay_field_select_string, logfileout='logs/msinfo.log')
    else:
        delay_state_select_string = ('STATE_ID in [%s'%delay_state_IDs[0])
        for state_ID in range(1,len(delay_state_IDs)):
            delay_state_select_string += (',%s')%delay_state_IDs[state_ID]
        delay_state_select_string += ']'
        subtable = tb.query(delay_state_select_string)
        delay_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        delay_scan_select_string = ','.join(["%s" % ii for ii in delay_scan_list])
        logprint ("Delay calibrator(s) scans are "+delay_scan_select_string, logfileout='logs/msinfo.log')
        delay_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        delay_field_select_string = ','.join(["%s" % ii for ii in delay_field_list])
        logprint ("Delay calibrator(s) are fields "+delay_field_select_string, logfileout='logs/msinfo.log')

    if (len(polarization_state_IDs) == 0):
        logprint ("No polarization calibration scans defined, no polarization calibration possible", logfileout='logs/msinfo.log')
        polarization_scan_select_string=''
        polarization_field_select_string=''
    else:
        logprint ("Warning: polarization calibration scans found, but polarization calibration not yet implemented", logfileout='logs/msinfo.log')
        polarization_state_select_string = ('STATE_ID in [%s'%polarization_state_IDs[0])
        for state_ID in range(1,len(polarization_state_IDs)):
            polarization_state_select_string += (',%s')%polarization_state_IDs[state_ID]
        polarization_state_select_string += ']'
        subtable = tb.query(polarization_state_select_string)
        polarization_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        polarization_scan_select_string = ','.join(["%s" % ii for ii in polarization_scan_list])
        logprint ("Polarization calibrator(s) scans are "+polarization_scan_select_string, logfileout='logs/msinfo.log')
        polarization_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        polarization_field_select_string = ','.join(["%s" % ii for ii in polarization_field_list])
        logprint ("Polarization calibrator(s) are fields "+polarization_field_select_string, logfileout='logs/msinfo.log')

    if (len(phase_state_IDs) == 0):
        QA2_msinfo='Fail'
        logprint("ERROR: No gain calibration scans found", logfileout='logs/msinfo.log')
        raise Exception("No gain calibration scans found")
    else:
        phase_state_select_string = ('STATE_ID in [%s'%phase_state_IDs[0])
        for state_ID in range(1,len(phase_state_IDs)):
            phase_state_select_string += (',%s')%phase_state_IDs[state_ID]
        phase_state_select_string += ']'
        subtable = tb.query(phase_state_select_string)
        phase_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        phase_scan_select_string = ','.join(["%s" % ii for ii in phase_scan_list])
        logprint ("Phase calibrator(s) scans are "+phase_scan_select_string, logfileout='logs/msinfo.log')
        phase_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        phase_field_select_string = ','.join(["%s" % ii for ii in phase_field_list])
        logprint ("Phase calibrator(s) are fields "+phase_field_select_string, logfileout='logs/msinfo.log')

    if (len(amp_state_IDs) == 0):
        logprint ("No amplitude calibration scans defined, will use phase calibrator", logfileout='logs/msinfo.log')
        amp_scan_select_string=phase_scan_select_string
        logprint ("Amplitude calibrator(s) scans are "+amp_scan_select_string, logfileout='logs/msinfo.log')
        amp_field_select_string=phase_scan_select_string
        logprint ("Amplitude calibrator(s) are fields "+amp_field_select_string, logfileout='logs/msinfo.log')
    else:
        amp_state_select_string = ('STATE_ID in [%s'%amp_state_IDs[0])
        for state_ID in range(1,len(amp_state_IDs)):
            amp_state_select_string += (',%s')%amp_state_IDs[state_ID]
        amp_state_select_string += ']'
        subtable = tb.query(amp_state_select_string)
        amp_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
        amp_scan_select_string = ','.join(["%s" % ii for ii in amp_scan_list])
        logprint ("Amplitude calibrator(s) scans are "+amp_scan_select_string, logfileout='logs/msinfo.log')
        amp_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
        amp_field_select_string = ','.join(["%s" % ii for ii in amp_field_list])
        logprint ("Amplitude calibrator(s) are fields "+amp_field_select_string, logfileout='logs/msinfo.log')

    # Find all calibrator scans and fields

    calibrator_state_select_string = ('STATE_ID in [%s'%calibrator_state_IDs[0])
    for state_ID in range(1,len(calibrator_state_IDs)):
        calibrator_state_select_string += (',%s')%calibrator_state_IDs[state_ID]

    calibrator_state_select_string += ']'
    subtable = tb.query(calibrator_state_select_string)
    calibrator_scan_list = list(np.unique(subtable.getcol('SCAN_NUMBER')))
    calibrator_scan_select_string = ','.join(["%s" % ii for ii in calibrator_scan_list])
    calibrator_field_list = list(np.unique(subtable.getcol('FIELD_ID')))
    calibrator_field_select_string = ','.join(["%s" % ii for ii in calibrator_field_list])

    tb.close()

if (((startdate >= 55918.80) and (startdate <= 55938.98)) or ((startdate >= 56253.6) and (startdate <= 56271.6))):
    logprint ("Weather station broken during this period, using 100%", logfileout='logs/msinfo.log')
    logprint ("seasonal model for calculating the zenith opacity", logfileout='logs/msinfo.log')
    tau=plotweather(vis=ms_active,seasonal_weight=1.0,doPlot=True)
else:
    tau=plotweather(vis=ms_active,seasonal_weight=0.5,doPlot=True)

#If there are any pointing state IDs, subtract 2 from the number of
#science spws -- needed for QA scores

if (len(pointing_state_IDs)>0):
    numSpws2 = numSpws - 2
else:
    numSpws2 = numSpws

logprint ("Zenith opacities based on weather data are: "+str(tau), logfileout='logs/msinfo.log')

#Prep string listing of correlations from dictionary created by method buildscans
#For now, only use the parallel hands.  Cross hands will be implemented later.
scandict = buildscans(ms_active)
corrstring_list = scandict['DataDescription'][0]['corrdesc']
removal_list = ['RL', 'LR', 'XY', 'YX']
corrstring_list = list(set(corrstring_list).difference(set(removal_list)))
corrstring = string.join(corrstring_list,',')
logprint ("Correlations shown in plotms will be "+corrstring, logfileout='logs/msinfo.log')

#Get number of antennas, store in numAntenna
tbLoc = casac.table()
tbLoc.open( '%s/ANTENNA' % ms_active)
nameAntenna = tbLoc.getcol( 'NAME' )
numAntenna = len(nameAntenna)
tbLoc.close()

minBL_for_cal=max(3,int(numAntenna/2.0))

#Determine if 3C84 was used as a bandpass or delay calibrator
positions = []

for ii in range(0,len(field_positions[0][0])):
    positions.append([field_positions[0][0][ii], field_positions[1][0][ii]])

fields_3C84 = find_3C84(positions)

cal3C84_d = False
cal3C84_bp = False
#uvrange3C84 = '0~1800klambda'
uvrange3C84 = ''

if fields_3C84 != []:
    for field_int in fields_3C84:
        if (str(field_int) in delay_field_select_string):
            cal3C84_d = True
            logprint("WARNING: 3C84 was observed as a delay calibrator, uvrange limits may be used", logfileout='logs/msinfo.log')

if fields_3C84 != []:
    for field_int in fields_3C84:
        if (str(field_int) in bandpass_field_select_string):
            cal3C84_bp = True
            logprint("WARNING: 3C84 was observed as a BP calibrator, uvrange limits may be used", logfileout='logs/msinfo.log')

# Identify bands/basebands/spws

tb.open(ms_active+'/SPECTRAL_WINDOW')
spw_names = tb.getcol('NAME')
tb.close()

# If the dataset is too old to have the bandname in it, assume that
# either there are 8 spws per baseband (and allow for one or two for
# pointing), or that this is a dataset with one spw per baseband

if (len(spw_names)>=8):
    critfrac=0.9/int(len(spw_names)/8.0)
else:
    critfrac=0.9/float(len(spw_names))

if '#' in spw_names[0]:
#
# i assume that if any of the spw_names have '#', they all do...
#
    bands_basebands_subbands = []
    for spw_name in spw_names:
        receiver_name, baseband, subband = spw_name.split('#')
        receiver_band = (receiver_name.split('_'))[1]
        bands_basebands_subbands.append([receiver_band, baseband, int(subband)])
    spws_info = [[bands_basebands_subbands[0][0], bands_basebands_subbands[0][1], [], []]]
    bands = [bands_basebands_subbands[0][0]]
    for ii in range(len(bands_basebands_subbands)):
        band,baseband,subband = bands_basebands_subbands[ii]
        found = -1
        for jj in range(len(spws_info)):
            oband,obaseband,osubband,ospw_list = spws_info[jj]
            if band==oband and baseband==obaseband:
                osubband.append(subband)
                ospw_list.append(ii)
                found = jj
                break
        if found >= 0:
            spws_info[found] = [oband,obaseband,osubband,ospw_list]
        else:
            spws_info.append([band,baseband,[subband],[ii]])
            bands.append(band)
    logprint("Bands/basebands/spws are:", logfileout='logs/msinfo.log')
    for spw_info in spws_info:
        spw_info_string = spw_info[0] + '   ' + spw_info[1] + '   [' + ','.join(["%d" % ii for ii in spw_info[2]]) + ']   [' + ','.join(["%d" % ii for ii in spw_info[3]]) + ']'
        logprint(spw_info_string, logfileout='logs/msinfo.log')
# Critical fraction of flagged solutions in delay cal to avoid an
# entire baseband being flagged on all antennas
    critfrac=0.9/float(len(spws_info))
elif ':' in spw_names[0]:
    logprint("old spw names with :", logfileout='logs/msinfo.log')
else:
    logprint("unknown spw names", logfileout='logs/msinfo.log')

# Check for missing scans

missingScans = 0
missingScanStr = ''

for i in range(max(scanNums)):
    if scanNums.count(i+1) == 1: pass
    else:
        logprint ("WARNING: Scan "+str(i+1)+" is not present", logfileout='logs/msinfo.log')
        missingScans += 1
        missingScanStr = missingScanStr+str(i+1)+', '

if (missingScans > 0):
    logprint ("WARNING: There were "+str(missingScans)+" missing scans in this MS", logfileout='logs/msinfo.log')
else:
    logprint ("No missing scans found.", logfileout='logs/msinfo.log')

# Make some plots and list ms contents

# Plot raw data

# (Needs to be done; Miriam has something in her script for this)

# Plot online flags

logprint ("Plotting online flags", logfileout='logs/msinfo.log')

syscommand='rm -rf onlineFlags.png'
os.system(syscommand)

default('flagcmd')
vis=ms_active
inpmode='xml'
action='plot'
plotfile='onlineFlags.png'
savepars=False
async=False
flagcmd()

logprint ("Online flags plot onlineFlags.png", logfileout='logs/msinfo.log')

# Plot antenna locations

logprint ("Plotting antenna positions", logfileout='logs/msinfo.log')

syscommand='rm -rf plotants.png'
os.system(syscommand)

default('plotants')
vis=ms_active
figfile='plotants.png'
async=False
plotants()

logprint ("Antenna positions plot plotants.png", logfileout='logs/msinfo.log')

# Make elevation vs. time plot

logprint ("Plotting elevation vs. time for all fields", logfileout='logs/msinfo.log')

syscommand='rm -rf el_vs_time.png'
os.system(syscommand)

#clearstat()

default('plotms')
vis=ms_active
xaxis='time'
yaxis='elevation'
selectdata=True
spw='*:31'
antenna='0&1;2&3'
correlation=corrstring
averagedata=False
transform=False
extendflag=False
iteraxis=''
customsymbol=False
coloraxis='field'
customflaggedsymbol=False
plotrange=[]
title='Elevation vs. time'
xlabel=''
ylabel=''
showmajorgrid=False
showminorgrid=False
plotfile='el_vs_time.png'
overwrite=True
showgui=False
async=False
plotms()


# attempt to fix missing plots in CASA 4.1.0 -- commented out for
# CASA 4.2.2 to see if it's still a problem
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
#
#spw='*:32'
#plotfile='el_vs_time.png'
#overwrite=True
#plotms()



logprint ("Elevation vs. time plot el_vs_time.png", logfileout='logs/msinfo.log')


logprint ("Finished EVLA_pipe_msinfo.py", logfileout='logs/msinfo.log')
logprint ("QA2 score: "+QA2_msinfo, logfileout='logs/msinfo.log')
time_list=runtiming('msinfo', 'end')

pipeline_save()

