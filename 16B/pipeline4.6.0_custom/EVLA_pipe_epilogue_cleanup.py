
'''
This script should be run once the calibration and flagging are final (by-hand
and via the pipeline).

Run this script to only keep the final cal tables, flag lists, and weblog.

The MS will also be split, keeping only the calibrated column and removing all
flagged data.

Final data products will be placed in a new folder, one level up, whose name
will be the current folder with "_calibrated" appended.
'''

import os
import shutil

try:
    ms_active
except NameError:
    raise Warning("Run EVLA_pipe_restore.py before this script.")

# Check if all non-essential pipeline products should be deleted
remove_nonessential = \
    True \
    if raw_input("Remove non-essential pipeline products (y/n): ") == "y" \
    else False

current_path = os.getcwd()

folder_name = current_path.split("/")[-1] + "_calibrated"
full_new_path = os.path.join(os.path.dirname(current_path), folder_name)

# Try making the new directory
if os.path.exists(full_new_path):
    raise IOError("The folder {} already exists. Remove and "
                  "re-run.".format(folder_name))

os.mkdir(full_new_path)

# Move final data products
shutil.move(os.path.join(current_path, "pipeline_shelf.restore"),
            full_new_path)
shutil.move(os.path.join(current_path, "final_caltables"), full_new_path)
shutil.move(os.path.join(current_path, "logs"), full_new_path)
shutil.move(os.path.join(current_path, "test_images"), full_new_path)
shutil.move(os.path.join(current_path, "weblog"), full_new_path)
shutil.move(os.path.join(current_path, ms_active + ".flagversions"),
            full_new_path)
# Some reductions have custom flagging scripts I've created. Ensure these are
# also kept.
shutil.move(os.path.join(current_path, "*.py"), full_new_path)

# Now split out the science fields, keep only the corrected column, and remove
# all flagged data.
default("split")
# For the continuum tracks, I want to keep everything for now to avoid
# potential issues when applying the polarization calibration
if "continuum" in ms_active:
    intent = ''
else:
    intent = "*OBSERVE*"

split(vis=ms_active,
      outputvis=os.path.join(full_new_path, ms_active[:-3] + "_calibrated.ms"),
      intent=intent, datacolumn='corrected', keepflags=False)

# Now nuke the old stuff
if remove_nonessential:
    os.chdir("..")
    os.system("rm -rf {}".format(current_path.split("/")[-1]))


# Now tar and zip
os.chdir(full_new_path)

os.system("tar -zcf weblog.tar.gz weblog")
os.system("rm -rf weblog")
os.system("tar -zcf logs.tar.gz logs")
os.system("rm -rf logs")
os.system("tar -zcf test_images.tar.gz test_images")
os.system("rm -rf test_images")
os.system("tar -zcf final_caltables.tar.gz final_caltables")
os.system("rm -rf final_caltables")
os.system("tar -zcf {0}.flagversions.tar.gz {0}.flagversions".format(ms_active))
os.system("rm -rf {0}.flagversions".format(ms_active))

os.mkdir("restoration_products")
shutil.move("weblog.tar.gz", "restoration_products")
shutil.move("logs.tar.gz", "restoration_products")
shutil.move("test_images.tar.gz", "restoration_products")
shutil.move("final_caltables.tar.gz", "restoration_products")
shutil.move("{}.flagversions.tar.gz".format(ms_active), "restoration_products")
shutil.move("pipeline_shelf.restore", "restoration_products")
shutil.move("*.py", "restoration_products")
