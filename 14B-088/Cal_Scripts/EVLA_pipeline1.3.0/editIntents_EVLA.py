
'''
VLA version of editIntents (https://casaguides.nrao.edu/index.php/EditIntents)
'''

import analysisUtils as au
import numpy as np


def editIntents(msName='', field='', scan='', newintents='', help=False,
                append=False):
    """
    Change the observation intents for a specified field.  Adapted from
    John Lightfoot's interactive function for the ALMA pipeline.
    For further help and examples, run editIntents(help=True) or
    see http://casaguides.nrao.edu/index.php?title=EditIntents
    - T. Hunter
    """
    validIntents = ['AMPLITUDE', 'ATMOSPHERE', 'BANDPASS', 'DELAY', 'FLUX',
                    'PHASE', 'SIDEBAND_RATIO', 'TARGET', 'WVR', 'POL_ANGLE',
                    'POL_LEAKAGE', 'CALIBRATE_AMPLI', 'SYS_CONFIG',
                    'CALIBRATE_ATMOSPHERE', 'CALIBRATE_BANDPASS',
                    'CALIBRATE_DELAY', 'CALIBRATE_FLUX', 'CALIBRATE_PHASE',
                    'CALIBRATE_SIDEBAND_RATIO', 'OBSERVE_TARGET',
                    'CALIBRATE_WVR', 'CALIBRATE_POL_ANGLE', 'CALIBRATE_POL_LEAKAGE']

    if help == True:
        print "Usage: editIntents(msName='', field='', scan='', newintents='', help=False, append=False)"
        print "valid intents = ", validIntents
        print "  Note: the first 9 are simply a shorthand way to specify the latter 9."
        print "  Enter as a list like above, or as single comma-delimited string."
        print "field is a required argument, while scan is optional."
        print "scan can be a single integer, integer list, or comma-delimited string list"
        print "Set append=True to add the specified intent to the existing intents."
        return
    if msName == '':
        print "You must specify a measurement set."
        return
    if field == '':
        print "You must specify a field ID or name."
        return

    mytb = au.createCasaTool(tbtool)
    mytb.open(msName + '/FIELD')
    fieldnames = mytb.getcol('NAME')
    print "Found fieldnames = ", fieldnames
    mytb.close()

    mytb.open(msName + '/STATE')
    intentcol = mytb.getcol('OBS_MODE')
    intentcol = intentcol
    mytb.close()

    mytb.open(msName, nomodify=False)
    naddedrows = 0
    if isinstance(newintents, list):
        desiredintents = ''
        for n in newintents:
            desiredintents += n
            if (n != newintents[-1]):
                desiredintents += ','
    else:
        desiredintents = newintents
    desiredintentsList = desiredintents.split(',')

    for intent in desiredintentsList:
        if intent not in validIntents:
            print "Invalid intent = %s.  Valid intents = %s" % (intent, validIntents)
            mytb.close()
            return
    foundField = False
    if (type(scan) != list):
        scan = str(scan)
    for id, name in enumerate(fieldnames):
      if name == field or id == field:
        foundField = True
        print 'FIELD_ID %s has name %s' % (id, name)
        if scan == '':
            s = mytb.query('FIELD_ID==%s' % id)
            print "Got %d rows in the ms matching field=%s" % (s.nrows(), id)
        else:
            if (type(scan) == str):
                scans = [int(x) for x in scan.split(',')]
            elif (type(scan) != list):
                scans = [scan]
            else:
                scans = scan
#            s = mytb.query('FIELD_ID==%s AND SCAN_NUMBER==%s' % (id, scan))
            print "Querying: 'FIELD_ID==%s AND SCAN_NUMBER in %s'" \
                % (id, str(scans))
            s = mytb.query('FIELD_ID==%s AND SCAN_NUMBER in %s'
                           % (id, str(scans)))
            print "Got %d rows in the ms matching field=%s and scan in %s" % \
                (s.nrows(), id, str(scans))
        if (s.nrows() == 0):
            mytb.close()
            print "Found zero rows for this field (and/or scan). Stopping."
            return
        state_ids = s.getcol('STATE_ID')
        # original code from J. Lightfoot, can probably be replaced
        # by the np.unique() above
        states = []
        for state_id in state_ids:
            if state_id not in states:
                states.append(state_id)
#        print "states = ", states
        for ij in range(len(states)):
            states[ij] = intentcol[states[ij]]
        print 'current intents are:'
        for state in states:
            print state

        if not append:
            states = []
        for desiredintent in desiredintentsList:
            if desiredintent.find('TARGET') >= 0:
                states.append('OBSERVE_TARGET#UNSPECIFIED')
            elif desiredintent.find('BANDPASS') >= 0:
                states.append('CALIBRATE_BANDPASS#UNSPECIFIED')
            elif desiredintent.find('PHASE') >= 0:
                states.append('CALIBRATE_PHASE#UNSPECIFIED')
            elif desiredintent.find('AMPLI') >= 0:
                states.append('CALIBRATE_AMPLI#UNSPECIFIED')
            elif desiredintent.find('FLUX') >= 0:
                states.append('CALIBRATE_FLUX#UNSPECIFIED')
            elif desiredintent.find('ATMOSPHERE') >= 0:
                states.append('CALIBRATE_ATMOSPHERE#UNSPECIFIED')
            elif desiredintent.find('WVR') >= 0:
                states.append('CALIBRATE_WVR#UNSPECIFIED')
            elif desiredintent.find('SIDEBAND_RATIO') >= 0:
                states.append('CALIBRATE_SIDEBAND_RATIO#UNSPECIFIED')
            elif desiredintent.find('DELAY') >= 0:
                states.append('CALIBRATE_DELAY#UNSPECIFIED')
            elif desiredintent.find('POL_ANGLE') >= 0:
                states.append('CALIBRATE_POL_ANGLE#UNSPECIFIED')
            elif desiredintent.find('POL_LEAKAGE') >= 0:
                states.append('CALIBRATE_POL_LEAKAGE#UNSPECIFIED')
            elif desiredintent.find('SYS_CONFIG') >= 0:
                states.append('SYSTEM_CONFIGURATION#UNSPECIFIED')
            else:
                print "Unrecognized intent = %s" % desiredintent
                continue
            print 'setting %s' % (states[-1])

        if states != []:
            state = reduce(lambda x, y: '%s,%s' % (x, y), states)
            if state not in intentcol:
                print 'adding intent to state table'
                intentcol = list(intentcol)
                intentcol.append(state)
                intentcol = np.array(intentcol)
                state_id = len(intentcol) - 1
                naddedrows += 1
                print 'state_id is', state_id
                state_ids[:] = state_id
            else:
                print 'intent already in state table'
                state_id = np.arange(len(intentcol))[intentcol == state]
                print 'state_id is', state_id
                if isinstance(state_id, list) or isinstance(state_id, np.ndarray):
                    # ms can have identical combinations of INTENT, so just
                    # pick the row for the first appearance - T. Hunter
                    state_ids[:] = state_id[0]
                else:
                    state_ids[:] = state_id
            s.putcol('STATE_ID', state_ids)
    if not foundField:
        print "Field not found"
        return
    mytb.close()

    print 'writing new STATE table'
    mytb.open(msName + '/STATE', nomodify=False)
    if naddedrows > 0:
        mytb.addrows(naddedrows)
    mytb.putcol('OBS_MODE', intentcol)
    mytb.close()
