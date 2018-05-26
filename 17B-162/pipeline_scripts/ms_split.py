
import sys
import os
from shutil import copyfile
from tasks import split, importasdm

'''
Identify the continuum and line SPWs and split into separate MSs and
directories.

Also searches for custom per-track flagging scripts in this repo.
'''

mySDM = sys.argv[-2]
parallel_run = True if sys.argv[-1] == "T" else False

ms_active = mySDM + ".ms"

print("Given inputs:")
print("SDM: {}".format(mySDM))
print("Make MMS: {}".format(parallel_run))

importasdm(asdm=mySDM, vis=ms_active, ocorr_mode='co',
           applyflags=True, savecmds=True, tbuff=1.5,
           outfile='{}.flagonline.txt'.format(mySDM),
           createmms=parallel_run)

parentdir = os.getcwd().split("/")[-1]

lines_folder = parentdir + '_speclines'
if not os.path.exists(lines_folder):
    os.mkdir(lines_folder)

# Copy the cont.dat file from the repo
source_path = os.path.expanduser("~/Dropbox/code_development/VLA_Lband/17B-162/pipeline_scripts/cont.dat")
copyfile(source_path, lines_folder + "/cont.dat")

# Check for a flag template and copy if in the repo
# Naming conventions is 17B-162_month_day_year_lines_flags.txt
flag_filename = "{}_lines_flags.txt".format(parentdir)
flag_path = os.path.expanduser("~/Dropbox/code_development/VLA_Lband/17B-162/pipeline_scripts/track_flagging")
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

# While it would be nice to remove the pol cal scans here, the pipeline
# will fail when running fluxboot because there is no other calibration
# field to transfer to! We can avoid this by just keeping all of the
# fields, even if they aren't used.
# field="0137+331=3C48,M33_2,M33_14,M33_6,M33_7_center,M33_12,M33_11,M33_8")

cont_folder = parentdir + '_continuum'
if not os.path.exists(cont_folder):
    os.mkdir(cont_folder)

# Check for a flag template and copy if in the repo
# Naming conventions is 17B-162_month_day_year_lines_flags.txt
flag_filename = "{}_continuum_flags.txt".format(parentdir)
flag_path = os.path.expanduser("~/Dropbox/code_development/VLA_Lband/17B-162/pipeline_scripts/track_flagging")
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
