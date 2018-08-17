
import sys
import os
from glob import glob
import shutil
import numpy as np

from tasks import plotms, flagdata, flagmanager

mySDM = sys.argv[-1]
myvis = mySDM if mySDM.endswith("ms") else mySDM + ".ms"

__rethrow_casa_exceptions = True
context = h_init()
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
context.set_state('ProjectSummary', 'proposal_code', '17B-162')
context.set_state('ProjectSummary', 'piname', 'Eric Koch')

try:
    hifv_importdata(ocorr_mode='co', nocopy=False, vis=[mySDM],
                    createmms='automatic', asis='Receiver CalAtmosphere',
                    overwrite=False)
    hifv_hanning(pipelinemode="automatic")

    # Completely flag the 3 SPWs that are RFI-dominated
    flagdata(vis=myvis, spw="0,1,4", flagbackup=False)

    # Flag regions channels that are continuously affected by RFI
    spw2 = "2:0~14;21~26;44~54;83~96"
    spw3 = "3:0~15"
    spw5 = "5:52~67"
    spw7 = "7:44~127"
    flag_str = ",".join([spw2, spw3, spw5, spw7])
    flagdata(vis=myvis, spw=flag_str, flagbackup=False)

    # Add extra quack for the first 3C48 scan and scans after longer slews
    # May need to be adjusted for this. Seem to need ~30 s more after online
    # flags indicate
    flagdata(vis=myvis, flagbackup=False, mode='quack', scan='2',
             quackinterval=270.0)
    flagdata(vis=myvis, flagbackup=False, mode='quack',
             scan='3,11,19,27,35,43,51,5967,75,83,91,99',
             quackinterval=14.0)
    flagdata(vis=myvis, flagbackup=False, mode='quack',
             scan='100',
             quackinterval=70.0)
    flagdata(vis=myvis, flagbackup=False, mode='quack',
             scan='100,101',
             quackinterval=100.0)

    flagmanager(vis=myvis, mode='save', versionname="known_RFI",
                comment="Removal of constant L-band RFI in 17B + Extra long-slew quacking.")

    # Online flags applied when importing ASDM
    hifv_flagdata(intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,*SIDEBAND_RATIO*, \
                  *UNKNOWN*, *SYSTEM_CONFIGURATION*, \
                  *UNSPECIFIED#UNSPECIFIED*',
                  flagbackup=False, scan=True, baseband=True, clip=True,
                  autocorr=True,
                  hm_tbuff='1.5int', template=True,
                  filetemplate="additional_flagging.txt",
                  online=False, tbuff=0.0,
                  fracspw=0.05, shadow=True, quack=True, edgespw=True)
    hifv_vlasetjy(fluxdensity=-1, scalebychan=True, reffreq='1GHz', spix=0)
    hifv_priorcals(tecmaps=False)
    hifv_testBPdcals(weakbp=False)
    hifv_flagbaddef(pipelinemode="automatic")
    hifv_checkflag(pipelinemode="automatic")
    hifv_semiFinalBPdcals(weakbp=False)
    hifv_checkflag(checkflagmode='semi')
    hifv_semiFinalBPdcals(weakbp=False)
    hifv_solint(pipelinemode="automatic")
    hifv_fluxboot(pipelinemode="automatic")
    hifv_finalcals(weakbp=False)
    hifv_applycals(flagdetailedsum=True, flagbackup=True, calwt=[True],
                   flagsum=True, gainmap=False)
    hifv_targetflag(intents='*CALIBRATE*,*TARGET*')
    hifv_statwt(pipelinemode="automatic")
    hifv_plotsummary(pipelinemode="automatic")
    hif_makeimlist(nchan=-1, calmaxpix=300, intent='PHASE,BANDPASS')
    hif_makeimages(tlimit=2.0, hm_negativethreshold=-999.0,
                   subcontms=False, hm_masking='none', masklimit=4,
                   maxncleans=1, hm_growiterations=-999, cleancontranges=False,
                   noise='1.0Jy', hm_minbeamfrac=-999.0, target_list={}, robust=-999.0,
                   parallel='automatic', weighting='briggs', hm_noisethreshold=-999.0,
                   hm_lownoisethreshold=-999.0, npixels=0, hm_sidelobethreshold=-999.0)

    # Make a folder of products for restoring the pipeline solution
    os.mkdir('products/')
    hifv_exportdata(products_dir='products/')

finally:
    h_save()

# Make a new directory for the imaging outputs
if not os.path.exists("image_outputs"):
    os.mkdir("image_outputs")

image_files = glob("oussid*")

for fil in image_files:
    shutil.move(fil, "image_outputs/")

# Now make a bunch of scan plots to make it easier to identify bad data
ms_active = mySDM

# SPWs to loop through
tb.open(os.path.join(ms_active, "SPECTRAL_WINDOW"))
spws = range(len(tb.getcol("NAME")))
nchans = tb.getcol('NUM_CHAN')
tb.close()

# Read the field names
tb.open(os.path.join(ms_active, "FIELD"))
names = tb.getcol('NAME')
numFields = tb.nrows()
tb.close()

# Intent names
tb.open(os.path.join(ms_active, 'STATE'))
intentcol = tb.getcol('OBS_MODE')
tb.close()

tb.open(ms_active)
scanNums = np.unique(tb.getcol('SCAN_NUMBER'))
field_scans = []
is_calibrator = np.empty_like(scanNums, dtype='bool')
is_all_flagged = np.empty((len(spws), len(scanNums)), dtype='bool')
for ii in range(numFields):
    subtable = tb.query('FIELD_ID==%s' % ii)
    field_scan = np.unique(subtable.getcol('SCAN_NUMBER'))
    field_scans.append(field_scan)

    # Is the intent for calibration?
    scan_intents = intentcol[np.unique(subtable.getcol("STATE_ID"))]
    is_calib = False
    for intent in scan_intents:
        if "CALIBRATE" in intent:
            is_calib = True
            break

    is_calibrator[field_scan - 1] = is_calib

    # Are any of the scans completely flagged?
    for spw in spws:
        for scan in field_scan:
            scantable = \
                tb.query("SCAN_NUMBER=={0} AND DATA_DESC_ID=={1}".format(scan,
                                                                         spw))
            if scantable.getcol("FLAG").all():
                is_all_flagged[spw, scan - 1] = True
            else:
                is_all_flagged[spw, scan - 1] = False

tb.close()

# Make folder for scan plots
scan_dir = "scan_plots"

if not os.path.exists(scan_dir):
    os.mkdir(scan_dir)

for spw_num in spws:
    print("On SPW {}".format(spw))

    spw_folder = os.path.join(scan_dir, "spw_{}".format(spw_num))
    if not os.path.exists(spw_folder):
        os.mkdir(spw_folder)
    else:
        # Make sure any old plots are removed first.
        os.system("rm {}/*.png".format(spw_folder))

    for ii in range(len(field_scans)):
        print("On field {}".format(names[ii]))
        for jj in field_scans[ii]:

            # Check if all of the data is flagged.
            if is_all_flagged[spw_num, jj - 1]:
                print("All data flagged in SPW {0} scan {1}".format(spw_num, jj))
                continue

            print("On scan {}".format(jj))

            # Amp vs. time
            default('plotms')
            vis = ms_active
            xaxis = 'time'
            yaxis = 'amp'
            ydatacolumn = 'corrected'
            selectdata = True
            field = names[ii]
            scan = str(jj)
            spw = str(spw_num)
            correlation = "RR,LL"
            averagedata = True
            avgbaseline = True
            transform = False
            extendflag = False
            plotrange = []
            title = 'Amp vs Time: Field {0} Scan {1}'.format(names[ii], jj)
            xlabel = ''
            ylabel = ''
            showmajorgrid = False
            showminorgrid = False
            plotfile = os.path.join(spw_folder, 'field_{0}_amp_scan_{1}.png'.format(names[ii], jj))
            overwrite = True
            showgui = False
            async = False
            plotms()

            # Amp vs. channel
            default('plotms')
            vis = ms_active
            xaxis = 'chan'
            yaxis = 'amp'
            ydatacolumn = 'corrected'
            selectdata = True
            field = names[ii]
            scan = str(jj)
            spw = str(spw_num)
            correlation = "RR,LL"
            averagedata = True
            avgbaseline = True
            transform = False
            extendflag = False
            plotrange = []
            title = 'Amp vs Chan: Field {0} Scan {1}'.format(names[ii], jj)
            xlabel = ''
            ylabel = ''
            showmajorgrid = False
            showminorgrid = False
            plotfile = os.path.join(spw_folder, 'field_{0}_amp_chan_scan_{1}.png'.format(names[ii], jj))
            overwrite = True
            showgui = False
            async = False
            plotms()

            # Plot amp vs uvdist
            default('plotms')
            vis = ms_active
            xaxis = 'uvdist'
            yaxis = 'amp'
            ydatacolumn = 'corrected'
            selectdata = True
            field = names[ii]
            scan = str(jj)
            spw = str(spw_num)
            # avgchannel = str(max(channels))
            # avgtime = '1e8s'
            correlation = "RR,LL"
            averagedata = True
            avgbaseline = False
            transform = False
            extendflag = False
            plotrange = []
            title = 'Amp vs UVDist: Field {0} Scan {1}'.format(names[ii], jj)
            xlabel = ''
            ylabel = ''
            showmajorgrid = False
            showminorgrid = False
            plotfile = os.path.join(spw_folder, 'field_{0}_amp_uvdist_scan_{1}.png'.format(names[ii], jj))
            overwrite = True
            showgui = False
            async = False
            plotms()

            if is_calibrator[jj - 1]:
                # Plot phase vs time
                default('plotms')
                vis = ms_active
                xaxis = 'time'
                yaxis = 'phase'
                ydatacolumn = 'corrected'
                selectdata = True
                field = names[ii]
                scan = str(jj)
                spw = str(spw_num)
                correlation = "RR,LL"
                averagedata = True
                avgbaseline = True
                transform = False
                extendflag = False
                plotrange = []
                title = 'Phase vs Time: Field {0} Scan {1}'.format(names[ii], jj)
                xlabel = ''
                ylabel = ''
                showmajorgrid = False
                showminorgrid = False
                plotfile = os.path.join(spw_folder, 'field_{0}_phase_scan_{1}.png'.format(names[ii], jj))
                overwrite = True
                showgui = False
                async = False
                plotms()

                # Plot phase vs channel
                default('plotms')
                vis = ms_active
                xaxis = 'chan'
                yaxis = 'phase'
                ydatacolumn = 'corrected'
                selectdata = True
                field = names[ii]
                scan = str(jj)
                spw = str(spw_num)
                correlation = "RR,LL"
                averagedata = True
                avgbaseline = True
                transform = False
                extendflag = False
                plotrange = []
                title = 'Phase vs Chan: Field {0} Scan {1}'.format(names[ii], jj)
                xlabel = ''
                ylabel = ''
                showmajorgrid = False
                showminorgrid = False
                plotfile = os.path.join(spw_folder, 'field_{0}_phase_chan_scan_{1}.png'.format(names[ii], jj))
                overwrite = True
                showgui = False
                async = False
                plotms()

                # Plot phase vs uvdist
                default('plotms')
                vis = ms_active
                xaxis = 'uvdist'
                yaxis = 'phase'
                ydatacolumn = 'corrected'
                selectdata = True
                field = names[ii]
                scan = str(jj)
                spw = str(spw_num)
                correlation = "RR,LL"
                avgchannel = "128"
                # avgtime = '1e8s'
                averagedata = True
                avgbaseline = False
                transform = False
                extendflag = False
                plotrange = []
                title = 'Phase vs UVDist: Field {0} Scan {1}'.format(names[ii], jj)
                xlabel = ''
                ylabel = ''
                showmajorgrid = False
                showminorgrid = False
                plotfile = os.path.join(spw_folder, 'field_{0}_phase_uvdist_scan_{1}.png'.format(names[ii], jj))
                overwrite = True
                showgui = False
                async = False
                plotms()