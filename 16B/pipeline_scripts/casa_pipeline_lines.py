import sys
import os
from glob import glob
import shutil
import numpy as np

from tasks import plotms

mySDM = sys.argv[-1]
myvis = mySDM if mySDM.endswith("ms") else mySDM + ".ms"

if not os.path.exists("cont.dat"):
    raise ValueError("The cont.dat file is not in the pipeline directory.")

__rethrow_casa_exceptions = True
context = h_init()
context.set_state('ProjectSummary', 'observatory',
                  'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
if "16B-236" in myvis:
    context.set_state('ProjectSummary', 'proposal_code', '16B-236')
elif "16B-242" in myvis:
    context.set_state('ProjectSummary', 'proposal_code', '16B-242')
else:
    raise ValueError("Could not determine project code from SDM name.")

context.set_state('ProjectSummary', 'piname', 'Eric Koch')

# A couple tracks in 16B-242 need to avoid the default refant
if "16B-242.sb32614458.eb32984330.57698.29168733796" in mySDM:
    refant_ignore = 'ea01'
elif "16B-242.sb32614458.eb32982501.57696.352311215276" in mySDM:
    refant_ignore = 'ea01,ea25'
elif "16B-242.sb32614458.eb32980749.57695.35484743056" in mySDM:
    refant_ignore = 'ea01,ea25'
else:
    refant_ignore = ""

try:
    hifv_importdata(ocorr_mode='co', nocopy=False, vis=[mySDM],
                    createmms='automatic', asis='Receiver CalAtmosphere',
                    overwrite=False)
# Hanning smoothing is turned off in the following step.
# In the case of extreme RFI, Hanning smoothing, however,
# may still be required.
# Avoid Hanning smoothing spectral lines
# hifv_hanning(pipelinemode="automatic")
# Online flags applied when importing ASDM
    hifv_flagdata(intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,*SIDEBAND_RATIO*, \
                          *UNKNOWN*, *SYSTEM_CONFIGURATION*, *UNSPECIFIED#UNSPECIFIED*',
                  flagbackup=False, scan=True, baseband=True, clip=True, autocorr=True,
                  hm_tbuff='1.5int', template=True,
                  filetemplate="additional_flagging.txt",
                  online=False, tbuff=0.0, fracspw=0.05,
                  shadow=True, quack=True, edgespw=True)

    # Beginning and end of each J0119+3210
    # Adjust for slew time at beginning
    flagdata(vis=myvis, mode='quack', field="J0119+3210", quackinterval=13,
             quackmode='beg', flagbackup=False)
    flagdata(vis=myvis, mode='quack', field="J0119+3210", quackinterval=4,
             quackmode='endb', flagbackup=False)

    # Pol cals
    flagdata(vis=myvis, mode='quack', field="J0319+4130", quackinterval=70,
             quackmode='beg', flagbackup=False)

    # It's usually a longer slew to 3C138. Go with 1:45 min
    flagdata(vis=myvis, mode='quack', field="3C138", quackinterval=105,
             quackmode='beg', flagbackup=False)

    # And finally the sources
    # Adjust for slew time at beginning
    scans_after_slew = "4,8,12,16,20,24,28,32,36,40,44,48,52,56,60"
    flagdata(vis=myvis, mode='quack', scan=scans_after_slew,
             quackinterval=8,
             quackmode='beg', flagbackup=False)
    scans_before_slew = "6,10,14,18,22,26,30,34,38,42,46,50,54,58,62"
    flagdata(vis=myvis, mode='quack', scan=scans_before_slew,
             quackinterval=6,
             quackmode='endb', flagbackup=False)

    flagmanager(vis=myvis, mode='save', versionname="extra_quacking",
                comment="Extra long-slew quacking.")

    hifv_vlasetjy(fluxdensity=-1, scalebychan=True, reffreq='1GHz', spix=0)
    hifv_priorcals(tecmaps=False)
    hifv_testBPdcals(weakbp=False, refantignore=refant_ignore)
    hifv_flagbaddef(pipelinemode="automatic")
    hifv_checkflag(pipelinemode="automatic")
    hifv_semiFinalBPdcals(weakbp=False, refantignore=refant_ignore)
    hifv_checkflag(checkflagmode='semi')
    hifv_semiFinalBPdcals(weakbp=False, refantignore=refant_ignore)
    hifv_solint(pipelinemode="automatic", refantignore=refant_ignore)
    hifv_fluxboot(pipelinemode="automatic", refantignore=refant_ignore)
    hifv_finalcals(weakbp=False, refantignore=refant_ignore)
    hifv_applycals(flagdetailedsum=True, flagbackup=True, calwt=[True],
                   flagsum=True, gainmap=False)
# Keep the following two steps in the script if cont.dat exists.
# Otherwise we recommend to comment out the next two tasks,
# or at least remove '*TARGET*' from the hifv_targetflag call
    hifv_targetflag(intents='*CALIBRATE*')
    hifv_statwt(pipelinemode="automatic")
    hifv_plotsummary(pipelinemode="automatic")
    hif_makeimlist(nchan=-1, calmaxpix=300, intent='PHASE,BANDPASS')
    hif_makeimages(tlimit=2.0, hm_negativethreshold=-999.0,
                   subcontms=False, hm_masking='none', masklimit=4,
                   maxncleans=1, hm_growiterations=-999, cleancontranges=False,
                   noise='1.0Jy', hm_minbeamfrac=-999.0, target_list={},
                   robust=-999.0, parallel='automatic', weighting='briggs',
                   hm_noisethreshold=-999.0, hm_lownoisethreshold=-999.0,
                   npixels=0, hm_sidelobethreshold=-999.0)

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

# Plot the bandpasses per SPW as well
bp_folder = "finalBPcal_plots"
if not os.path.exists(bp_folder):
    os.mkdir(bp_folder)

tb.open(ms_active + "/SPECTRAL_WINDOW")
nspws = tb.getcol("NAME").shape[0]
tb.close()

for ii in range(nspws):
    filename = 'finalBPcal_amp_spw_' + str(ii) + '.png'
    syscommand = 'rm -rf ' + filename
    os.system(syscommand)

    default('plotcal')
    caltable = mySDM + '.finalBPcal.b'
    xaxis = 'freq'
    yaxis = 'amp'
    poln = ''
    field = ''
    antenna = ''
    spw = str(ii)
    timerange = ''
    subplot = 111
    overplot = False
    clearpanel = 'Auto'
    iteration = ''
    showflags = False
    plotsymbol = 'o'
    plotcolor = 'blue'
    markersize = 5.0
    fontsize = 10.0
    showgui = False
    figfile = os.path.join(bp_folder, filename)
    async = False
    plotcal()

for ii in range(nspws):
    filename = 'finalBPcal_phase_spw_' + str(ii) + '.png'
    syscommand = 'rm -rf ' + filename
    os.system(syscommand)

    antPlot = str(ii * 3) + '~' + str(ii * 3 + 2)

    default('plotcal')
    caltable = mySDM + '.finalBPcal.b'
    xaxis = 'freq'
    yaxis = 'phase'
    poln = ''
    field = ''
    antenna = ''
    spw = str(ii)
    timerange = ''
    subplot = 111
    overplot = False
    clearpanel = 'Auto'
    iteration = ''
    # plotrange=[0,0,-phaseplotmax,phaseplotmax]
    showflags = False
    plotsymbol = 'o'
    plotcolor = 'blue'
    markersize = 5.0
    fontsize = 10.0
    showgui = False
    figfile = os.path.join(bp_folder, filename)
    async = False
    plotcal()

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

    # Plotting the HI spw (0) takes so so long.
    # Make some simplifications to save time
    if spw_num == 0:
        avg_chan = "4"
    else:
        avg_chan = "1"

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
                print("All data flagged in SPW {0} scan {1}"
                      .format(spw_num, jj))
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
            avgchannel = str(avg_chan)
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
            plotfile = os.path.join(
                spw_folder, 'field_{0}_amp_scan_{1}.png'.format(names[ii], jj))
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
            avgchannel = str(avg_chan)
            avgtime = "1e8"
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
            plotfile = os.path.join(
                spw_folder, 'field_{0}_amp_chan_scan_{1}.png'.format(names[ii], jj))
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
            avgchannel = str(4096)
            avgtime = '1e8'
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
            plotfile = os.path.join(
                spw_folder, 'field_{0}_amp_uvdist_scan_{1}.png'.format(names[ii], jj))
            overwrite = True
            showgui = False
            async = False
            plotms()

            # Skip the phase plots for the HI SPW (0)
            if is_calibrator[jj - 1]:
                # Plot phase vs time
                default('plotms')
                vis = ms_active
                xaxis = 'time'
                yaxis = 'phase'
                ydatacolumn = 'corrected'
                avgchannel = str(4096)
                selectdata = True
                field = names[ii]
                scan = str(jj)
                spw = str(spw_num)
                correlation = "RR,LL"
                averagedata = True
                avgbaseline = False
                transform = False
                extendflag = False
                plotrange = []
                title = 'Phase vs Time: Field {0} Scan {1}'.format(
                    names[ii], jj)
                xlabel = ''
                ylabel = ''
                showmajorgrid = False
                showminorgrid = False
                plotfile = os.path.join(
                    spw_folder, 'field_{0}_phase_scan_{1}.png'.format(names[ii], jj))
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
                avgchannel = str(avg_chan)
                correlation = "RR,LL"
                averagedata = True
                avgbaseline = False
                transform = False
                extendflag = False
                plotrange = []
                title = 'Phase vs Chan: Field {0} Scan {1}'.format(
                    names[ii], jj)
                xlabel = ''
                ylabel = ''
                showmajorgrid = False
                showminorgrid = False
                plotfile = os.path.join(
                    spw_folder, 'field_{0}_phase_chan_scan_{1}.png'.format(names[ii], jj))
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
                avgchannel = "4096"
                avgtime = '1e8'
                averagedata = True
                avgbaseline = False
                transform = False
                extendflag = False
                plotrange = []
                title = 'Phase vs UVDist: Field {0} Scan {1}'.format(
                    names[ii], jj)
                xlabel = ''
                ylabel = ''
                showmajorgrid = False
                showminorgrid = False
                plotfile = os.path.join(
                    spw_folder, 'field_{0}_phase_uvdist_scan_{1}.png'.format(names[ii], jj))
                overwrite = True
                showgui = False
                async = False
                plotms()
