
'''
Given a SPW, split it out from the set of tracks. Assumes each track has the
same setup.
'''

import glob
import sys
import os

from casa_tools import mysplit

# Set the project code to
proj_code = "14B-088"

#  SPW setup for 14B-088
spw_dict = {1: ["HI", "1.420405752GHz"],
            4: ["OH1612", "1.612231GHz"],
            6: ["OH1665", "1.6654018GHz"],
            7: ["OH1667", "1.667359GHz"],
            8: ["OH1720", "1.72053GHz"],
            0: ["H172alp", "1.28117GHz"],
            10: ["H152alp", "1.85425GHz"],
            9: ["H153alp", "1.81825GHz"],
            2: ["H166alp", "1.42473GHz"],
            5: ["H158alp", "1.65154GHz"],
            3: ["H164alp", "1.47734GHz"]}

# SPW
spw = int(sys.argv[-4])
intent = sys.argv[-3]
input_dir = sys.argv[-2]
output_dir = sys.argv[-1]

if intent == "None":
    intent = ""
elif intent.lower() == "target":
    intent = "[OBSERVE_TARGET#UNSPECIFIED]"
else:
    raise Warning("Only supporting intents of 'None' (ie. all) or 'target'")

spw_outdir = os.path.join(output_dir, spw_dict[spw][0])

# Check if the output directory exists
if not os.path.isdir(spw_outdir):
    os.mkdir(spw_outdir)


# Grab the folders containing the tracks. I've arranged them with the
# observations date, which will be used in the name of the split MSs.
tracks = glob.glob(os.path.join(input_dir, proj_code+"*lines"))

print("Tracks found: ")
for track in tracks:
    print(track)

print("Splitting out SPW "+str(spw)+" corresponding to "+spw_dict[spw][0])

for i, track in enumerate(tracks):
    print("Splitting from "+track.rstrip("/").split("/")[-1])
    print(str(i+1)+" out of "+str(len(tracks)))

    # Includes 14B-088 prefix
    obs_date = \
        track.rstrip("/").split("/")[-1].rstrip("_lines").lstrip(proj_code)

    split_name = proj_code+"_"+obs_date+"_"+spw_dict[spw][0]+"_" + \
        spw_dict[spw][1]+".ms"

    full_name = os.path.join(spw_outdir, split_name)

    # Now we have to get the actual ms name out of the pipeline output
    # My folder structures are a bit different due to my custom pipeline
    track_contents = os.listdir(track)

    # The tracks are too big to rewrite right now, so some still have
    # the 'extra_flagging' directory with pipeline re-runs. Check if we
    # should be looking for that.
    if "extra_flagging" in track_contents:
        track = os.path.join(track, "extra_flagging")
        track_contents = os.listdir(track)

    # Pick out the ms that starts out with the project code.
    ms_name = None
    posn = 0
    while posn < len(track_contents):
        if track_contents[posn].startswith(proj_code):
            if track_contents[posn].endswith(".ms"):
                ms_name = track_contents[posn]
                break
        posn += 1
    else:
        raise Warning("Could not find the ms in "+track)

    inputvis = os.path.join(track, ms_name)

    mysplit(vis=inputvis, outputvis=full_name, spw=str(spw),
            datacolumn="CORRECTED", keepflags=False, intent=intent)
