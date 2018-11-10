'''
The track taken on 11/05/16 was split in two due to a malfunction.
Concatenate the two parts together before splitting into line and
continuum MS.
'''


import sys
import os
from shutil import copyfile
import numpy as np
from tasks import split, importasdm, concat

mySDM = "16B-242.sb32614458.eb32984320.57697.291263148145"
mySDM2 = "16B-242_sb32614458_5.57697.30203707176"

parallel_run = True if sys.argv[-2] == "T" else False

# Split out the lines, continuum or both
which_split = sys.argv[-1]
if which_split == 'lines':
    split_lines = True
    split_cont = False
elif which_split == 'cont':
    split_cont = True
    split_lines = False
elif which_split == 'all':
    split_cont = True
    split_lines = True
else:
    raise ValueError("Input for splitting should be lines, "
                     "cont or all. Given {}".format(which_split))


ms_active = mySDM + ".ms"
ms_active2 = mySDM2 + ".ms"

print("Given inputs:")
print("SDM: {}".format(mySDM))
print("Make MMS: {}".format(parallel_run))

project_code = "16B-242"

if not os.path.exists(ms_active) and not os.path.exists(ms_active2):
    importasdm(asdm=mySDM, vis=ms_active, ocorr_mode='co',
               applyflags=True, savecmds=True, tbuff=1.5,
               outfile='{}.flagonline.txt'.format(mySDM),
               createmms=parallel_run)

    importasdm(asdm=mySDM2, vis=ms_active2, ocorr_mode='co',
               applyflags=True, savecmds=True, tbuff=1.5,
               outfile='{}.flagonline.txt'.format(mySDM2),
               createmms=parallel_run)

    # Concat these together
    concat(vis=[ms_active, ms_active2], outputvis=ms_active2,
           timesort=True)

    # But the pipeline will treat these as individual observations
    # We don't want that. So we need to remove information from the
    # ObsID
    tb.open(ms_active, nomodify=False)
    obsid = tb.getcol('OBSERVATION_ID')
    tb.putcol('OBSERVATION_ID', np.zeros_like(obs_ids))
    tb.close()

    # We also need to alter the OBSERVATION table
    tb.open(ms_active + "/OBSERVATION", nomodify=False)

    time_range = tb.getcol("TIME_RANGE")
    new_time_range = np.array([time_range.min(),
                               time_range.max()]).reshape((2, 1))
    new_log = tb.getcol("LOG")[0][:1].reshape((1, 1))

    # Delete the second column
    tb.removerows(1)

    tb.putcol("TIME_RANGE", new_time_range)

    tb.close()


else:
    print("MS already exists. Skipping importasdm")

parentdir = os.getcwd().split("/")[-1]

if split_lines:

    lines_folder = parentdir + '_speclines'
    if not os.path.exists(lines_folder):
        os.mkdir(lines_folder)

    project_path = os.path.expanduser("~/ownCloud/code_development/VLA_Lband/16B/pipeline_scripts")

    # Copy the cont.dat file from the repo
    if project_code == "16B-236":
        source_path = os.path.join(project_path, "cont_236.dat")
    else:
        source_path = os.path.join(project_path, "cont_242.dat")
    copyfile(source_path, lines_folder + "/cont.dat")

    # Check for a flag template and copy if in the repo
    # Naming conventions is 16B-236_month_day_year_lines_flags.txt
    flag_filename = "{}_lines_flags.txt".format(parentdir)
    flag_path = os.path.expanduser("~/ownCloud/code_development/VLA_Lband/16B/{}/track_flagging".format(project_code))
    full_flag_filename = os.path.join(flag_path, flag_filename)

    if os.path.exists(full_flag_filename):
        copyfile(full_flag_filename,
                 os.path.join(lines_folder, "additional_flagging.txt"))
    else:
        print("No additional flagging script found in the VLA_Lband repo"
              " for lines.")

    split(vis=ms_active,
          outputvis=lines_folder + "/" + mySDM + ".speclines.ms",
          spw="8~17", datacolumn='DATA', field="")

if split_cont:

    cont_folder = parentdir + '_continuum'
    if not os.path.exists(cont_folder):
        os.mkdir(cont_folder)

    # Check for a flag template and copy if in the repo
    # Naming conventions is 17B-162_month_day_year_lines_flags.txt
    flag_filename = "{}_continuum_flags.txt".format(parentdir)
    flag_path = os.path.expanduser("~/ownCloud/code_development/VLA_Lband/16B/{}/track_flagging".format(project_code))
    full_flag_filename = os.path.join(flag_path, flag_filename)
    if os.path.exists(full_flag_filename):
        copyfile(full_flag_filename,
                 os.path.join(cont_folder, "additional_flagging.txt"))
    else:
        print("No additional flagging script found in the VLA_Lband repo"
              " for continuum.")

    split(vis=ms_active,
          outputvis=cont_folder + "/" + mySDM + ".continuum.ms",
          spw="0~7", datacolumn='DATA',
          field="")
