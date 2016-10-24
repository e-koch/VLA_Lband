
'''
Execute with:
casa -c find_flag_badscans.py True True
within a directory with the pipeline products and MS.
'''

import os
import matplotlib.pyplot as p
import matplotlib.image as im
import json
from datetime import datetime
from glob import glob
import numpy as np

# Enable interactive
p.ion()


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


find_bad_scans = \
    True if raw_input("Bad scan finding? (True/False): ") == "True" else False
plot_bad_scans = \
    True if raw_input("Bad scan plotting? (True/False): ") == "True" else False

try:
    ms_active
except NameError:
    ms_active = raw_input("Provide the MS name: ")

out_file = os.path.join(ms_active + "_badscans.json")

bad_scan_dict = {}

if find_bad_scans:
    # Check if the scan plots directory exists
    plot_direc = ms_active[:-3] + "_scan_plots"
    if not os.path.exists(plot_direc):
        raise Exception("Create the scan plots directory first!")

    spw_folders = listdir_fullpath(plot_direc)

    print("Mark 'T' or 't' for plotting later, nothing to skip, and "
          "'A' or 'a' to flag the whole scan.")

    for spw_folder in spw_folders:

        bad_scan_dict[spw_folder.split("/")[-1]] = {}

        scan_plots = listdir_fullpath(spw_folder)

        # Keep track of scans that should be completely flagged.
        all_bad_scan = ""

        bad_scans_amp = ""

        amp_scan_plots = \
            np.array([plot for plot in scan_plots if "amp" in plot])

        # Now sort by the scan number
        amp_scan_nums = \
            np.array([int(plot.split("/")[-1].split("_")[-1].rstrip(".png"))
                      for plot in amp_scan_plots])
        amp_scan_plots = amp_scan_plots[np.argsort(amp_scan_nums)]

        print("Plotting amp vs time plots.")
        for plot in amp_scan_plots:
            num = plot.split("/")[-1].split("_")[-1].rstrip(".png")
            img = im.imread(plot)
            p.imshow(img, origin='upper')
            p.axis("off")
            inp = raw_input("Scan %s T(t)/''(f)/A(a): " % (num))
            bad = True if inp == "T" or inp == "t" else False
            if inp == "T" or inp == 't':
                bad_scans_amp += num + ","
            elif inp == "A" or inp == 'a':
                all_bad_scan += num + ","
            else:
                pass
            p.clf()

        bad_scans_amp = bad_scans_amp.rstrip(",")

        bad_scan_dict[spw_folder.split("/")[-1]]["Amp"] = bad_scans_amp

        phase_scan_plots = \
            np.array([plot for plot in scan_plots if "phase" in plot])

        if len(phase_scan_plots) != 0:

            # Now sort by the scan number
            phase_scan_nums = \
                np.array([int(plot.split("/")[-1].split("_")[-1].rstrip(".png"))
                          for plot in phase_scan_plots])
            phase_scan_plots = phase_scan_plots[np.argsort(phase_scan_nums)]

            print("Plotting phase vs time plots.")

            bad_scans_phase = ""

            for plot in phase_scan_plots:
                num = plot.split("/")[-1].split("_")[-1].rstrip(".png")
                img = im.imread(plot)
                p.imshow(img, origin='upper')
                p.axis("off")
                inp = raw_input("Scan %s T(t)/''(f)/A(a): " % (num))
                if inp == "T" or inp == 't':
                    bad_scans_phase += num + ","
                elif inp == "A" or inp == 'a':
                    all_bad_scan += num + ","
                else:
                    pass
                p.clf()

            bad_scans_phase = bad_scans_phase.rstrip(",")

            bad_scan_dict[spw_folder.split("/")[-1]]["Phase"] = bad_scans_phase

        all_bad_scan = all_bad_scan.rstrip(",")

        bad_scan_dict[spw_folder.split("/")[-1]]["Flag"] = all_bad_scan

    # Now save the dictionary via json
    with open(out_file, "a") as f:
        json.dump(bad_scan_dict, f)

# If enabled, now plot those bad scans to be flagged
# (Obviously this should be run in CASA to do this).
if plot_bad_scans:

    try:
        from tasks import plotms, flagdata, flagmanager
    except ImportError:
        raise ImportError("Run in CASA to plot!")

    # Check the bad scans file exists
    if not os.path.exists(out_file):
        raise Exception("No bad scans file found. "
                        "Enable find_bad_scans first.")

    # Load the bad scan file
    with open(out_file, "r") as f:
        bad_scan_dict = json.load(f)

    for spw in bad_scan_dict:

        spw_num = str(spw[-1])

        print("Plotting SPW {}".format(spw_num))

        bad_scans = bad_scan_dict[spw]

        # Check for scans that need to be entirely flagged.
        if len(bad_scans["Flag"]) != 0:
            print("Found the following scans to flag:"
                  " {}".format(bad_scans["Flag"]))

            default("flagdata")
            flagdata(vis=ms_active, spw=spw_num,
                     scan=bad_scans["Flag"], flagbackup=False)

        if len(bad_scans["Amp"]) != 0:
            print("The bad scans are {}".format(bad_scans["Amp"]))

            default('plotms')
            vis = ms_active
            xaxis = 'time'
            yaxis = 'amp'
            ydatacolumn = 'corrected'
            selectdata = True
            field = ''
            scan = str(bad_scans["Amp"])
            spw = spw_num
            correlation = "RR,LL"
            iteraxis = "scan"
            averagedata = True
            avgbaseline = True
            transform = False
            extendflag = False
            plotrange = []
            xlabel = ''
            ylabel = ''
            showmajorgrid = False
            showminorgrid = False
            overwrite = False
            showgui = True
            async = False
            plotms()
        else:
            print("No bad amp vs. time scans.")

        if "Phase" in bad_scans.keys():
            if len(bad_scans["Phase"]) != 0:
                print("The bad scans are {}".format(bad_scans["Phase"]))

                default('plotms')
                vis = ms_active
                xaxis = 'time'
                yaxis = 'phase'
                ydatacolumn = 'corrected'
                selectdata = True
                field = ''
                scan = str(bad_scans["Phase"])
                spw = spw_num
                correlation = "RR,LL"
                iteraxis = "scan"
                averagedata = True
                avgbaseline = True
                transform = False
                extendflag = False
                plotrange = []
                xlabel = ''
                ylabel = ''
                showmajorgrid = False
                showminorgrid = False
                overwrite = False
                showgui = True
                async = False
                plotms()

            else:
                print("No bad phase vs. time scans.")

        raw_input("Continue?")

    # Get the existing flag version names.
    flag_folder = "{}.flagversions".format(ms_active)
    tstamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    if not os.path.exists(flag_folder):
        print("No flag versions exist. Using default flag name.")
        versionname = "badscan_flagging_1_{}".format(tstamp)
    else:
        flag_versions = \
            glob(os.path.join(flag_folder, "flags.badscan_flagging_*"))
        if len(flag_versions) == 0:
            versionname = "badscan_flagging_1_{}".format(tstamp)
        else:
            num = len(flag_versions) + 1
            versionname = "badscan_flagging_{0}_{1}".format(num, tstamp)

    # Save this new version of the flags
    flagmanager(vis=ms_active, mode='save', versionname=versionname)
