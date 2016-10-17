
import sys
import os
import numpy as np

'''
Easier searching for good RFI flagging values
'''


def if_empty_return_old(string, old_val):
    if string == "":
        return old_val
    else:
        return float(string)


try:
    ms_name = sys.argv[1]
    apply_flagging = True if sys.argv[2] == "True" else False
    extend_pol = True if sys.argv[3] == "True" else False
except IndexError:
    ms_name = raw_input("Input vis? : ")
    apply_flagging = \
        True if raw_input("Apply the flagging? : ") == "True" else False
    extend_pol = \
        True if raw_input("Extend across pols? : ") == "True" else False

# Just want the number of SPWs
tb.open(os.path.join(ms_name, "SPECTRAL_WINDOW"))
nchans = tb.getcol('NUM_CHAN')
tb.close()

spws = range(len(nchans))

params_used = np.empty((len(nchans), 4))

default('flagdata')

for spw in spws:
    print("On spw "+str(spw)+" of "+str(len(nchans)))
    freqdevscale = 4.0
    timedevscale = 4.0
    growtime = 99.0
    growfreq = 99.0
    print("Starting at freqdevscale = %s and timedevscale = %s" %
          (freqdevscale, timedevscale))
    print("Starting at growfreq = %s and growtime = %s" % (growfreq, growtime))
    while True:
        flagdata(vis=ms_name, mode='rflag', field='*3C48*',
                 spw=str(spw), datacolumn='corrected',
                 action='calculate', display='both',
                 freqdevscale=freqdevscale, timedevscale=timedevscale,
                 growtime=growtime, growfreq=growfreq,
                 flagbackup=False)

        adjust = True if raw_input("New thresholds? : ") == "True" else False

        if adjust:
            print("Current freqdevscale and timedevscale: %s %s" %
                  (freqdevscale, timedevscale))
            freqdevscale = \
                if_empty_return_old(raw_input("New freqdevscale : "),
                                    freqdevscale)
            timedevscale = \
                if_empty_return_old(raw_input("New timedevscale : "),
                                    timedevscale)
            growfreq = if_empty_return_old(raw_input("New growfreq : "),
                                           growfreq)
            growtime = if_empty_return_old(raw_input("New growtime : "),
                                           growtime)
        else:
            break

    # Now apply the flagging
    if apply_flagging:

        apply_rfi = True if raw_input("Apply RFI flagging? : ") == "True" else False

        if apply_rfi:

            if extend_pol:
                flagdata(vis=ms_name, mode='rflag', spw=str(spw),
                         freqdevscale=freqdevscale, timedevscale=timedevscale,
                         growfreq=growfreq, growtime=growtime,
                         datacolumn='corrected',
                         action='apply', display='report', flagbackup=False)
                flagdata(vis=ms_name, mode='extend', extendpols=True,
                         spw=str(spw), action='apply', display='report',
                         datacolumn='corrected',
                         flagbackup=False)
            else:
                flagdata(vis=ms_name, spw=str(spw), mode='rflag',
                         freqdevscale=freqdevscale, timedevscale=timedevscale,
                         growfreq=growfreq, growtime=growtime,
                         datacolumn='corrected',
                         extendpols=False, action='apply', display='report',
                         flagbackup=False)

            params_used[spw, :] = \
                [freqdevscale, timedevscale, growfreq, growtime]

        else:
            params_used[spw, :] = \
                [np.NaN, np.NaN, np.NaN, np.NaN]

        obliterate = \
            True if raw_input("Flag whole SPW? : ") == "True" else False

        if obliterate:
            flagdata(vis=ms_name, spw=str(spw), mode='apply')
            print("Flagged all of SPW " + str(spw))
            params_used[spw, :] = [0.0, 0.0, 0.0, 0.0]

    else:
        print("Applying not enabled.")

if apply_flagging:
    np.savetxt(ms_name[:-3]+"_RFI_params.txt", params_used,
               header="freqdevscale, timedevscale, growfreq, growtime")

    flagmanager(vis=ms_name, mode='save',
                versionname='after_manual_rflagging_1')

