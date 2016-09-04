
'''
A few of the scan intents are wrong. This script keeps track of the changes
that are needed, based on the default SDM names.
'''

from CASA_functions import editIntents


# If run in a pipeline environment, try to read the MS name. Otherwise input
# needed

try:
    ms_active
except NameError:
    ms_active = raw_input("Give MS name to change intents: ")
    ms_active = ms_active.rstrip("/")

if "12A-403.sb12393368.eb13591366.56218.123318067126" in ms_active:

    # Setup
    editIntents(msName=ms_active, field='3C48', scan='0',
                newintents=['SYS_CONFIG'])

    # Set the BP, FLUX, DELAY cal
    editIntents(msName=ms_active, field='3C48', scan='2,3',
                newintents=['BANDPASS', 'DELAY', 'FLUX'])

    # Set the GAIN cal
    editIntents(msName=ms_active, field='J0119+3210', scan='',
                newintents=['AMPLITUDE', 'PHASE'])

    # Set the POL ANGLE and LEAKAGE
    # Include TARGET so the pipeline applies calibration to pol cal (otherwise
    # it isn't used by default).
    editIntents(msName=ms_active, field='3C84', scan='',
                newintents=['POL_ANGLE', 'POL_LEAKAGE', "TARGET"])
else:
    raise NameError("Given MS name ({}) does not match any of the SB names"
                    " from this project.")
