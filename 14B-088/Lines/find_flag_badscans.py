
# First create the file structure and splits with split_spw_from_all_tracks.py
# then make the scan plots using plot_scans.py

# Execute with:
# casa -c find_flag_badscans.py directory_with_MSs True True

import sys
import os
import matplotlib.pyplot as p
import matplotlib.image as im
import json

# Enable interactive
p.ion()


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

line_direc = sys.argv[-3]
find_bad_scans = True if sys.argv[-2] == "True" else False
plot_bad_scans = True if sys.argv[-1] == "True" else False

line_name = line_direc.rstrip("/").split("/")[-1]

out_file = os.path.join(line_direc, line_name + "_badscans.txt")

bad_scan_dict = {}
ms_names = []

if find_bad_scans:
    for f in listdir_fullpath(line_direc):
        if f.endswith(".ms"):
            print("Found MS!: " + f.split("/")[-1])
            # Check if the scan plots directory exists
            plot_direc = os.path.join(line_direc, "scan_plots",
                                      f.rstrip(".ms")+"_scan_plots")
            if not os.path.exists(plot_direc):
                raise Exception("Create the scan plots directory first!")

            ms_names.append(f)
            scan_plots = listdir_fullpath(plot_direc)

            bad_scans = ""

            for plot in scan_plots:
                num = plot.split("/")[-1].split("_")[-1].rstrip(".png")
                img = im.imread(plot)
                p.imshow(img, origin='upper')
                p.axis("off")
                inp = raw_input("Scan %s : " % (num))
                bad = True if inp == "T" or inp == "t" else False
                if bad:
                    bad_scans += num + ","
                p.clf()

            bad_scan_dict[f.split("/")[-1]] = bad_scans.rstrip(",")

    # Now save the dictionary via json
    with open(out_file, "a") as f:
        json.dump(bad_scan_dict, f)
else:
    # Check the bad scans file exists
    if not os.path.exists(out_file):
        raise Exception("No bad scans file found. "
                        "Enable find_bad_scans first.")

    # Load the bad scan file
    with open(out_file, "r") as f:
        bad_scan_dict = json.load(f)

    ms_names = []

    for f in listdir_fullpath(line_direc):
        if f.endswith(".ms"):
            if not f.split("/")[-1] in [str(key) for key in bad_scan_dict.keys()]:
                raise KeyError("Bad scan file has no MS name: " +
                               f.split("/")[-1])
            else:
                ms_names.append(str(f))


# If enabled, now plot those bad scans to be flagged
# (Obviously this should be run in CASA to do this).
if plot_bad_scans:
    try:
        from tasks import plotms
    except ImportError:
        raise ImportError("Run in CASA to plot!")

    for ms in ms_names:

        print("Plotting " + ms.split("/")[-1])

        default('plotms')
        vis = ms
        xaxis = 'time'
        yaxis = 'amp'
        ydatacolumn = 'corrected'
        selectdata = True
        field = ''
        scan = str(bad_scan_dict[ms.split("/")[-1]]).rstrip(",")
        spw = ''
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

        raw_input("Continue?")
