
'''
Backup an MS's flag backup and the final calibration tables.
'''

import sys
import os
from distutils.dir_util import copy_tree
from distutils.file_util import copy_file
import shelve
import shutil
from glob import glob

from tasks import flagmanager


def pipeline_restore(shelf_filename='pipeline_shelf.restore'):
    '''Restore the state of the pipeline from shelf file
    '''
    if os.path.exists(shelf_filename):
        try:
            pipe_shelf = shelve.open(shelf_filename)
        except Exception, e:
            logprint ("Restore point does not exist: "+str(e))

        for key in pipe_shelf:
            try:
                globals()[key] = pipe_shelf[key]
                key_status = True
            except:
                key_status = False


ms_path = sys.argv[-2]
out_path = sys.argv[-1]

# Restore the pipeline variables
pipeline_restore(shelf_filename=os.path.join(ms_path,
                                             "pipeline_shelf.restore"))

ms_active = os.path.join(ms_path, ms_active)

if not os.path.exists(out_path):
    os.mkdir(out_path)

# Save the final flag versions.
# default('flagmanager')
flagmanager(vis=ms_active, mode='save',
            versionname='finalflags_wmanualflagging',
            comment='Final flags after all calibration and flagging',
            merge='replace')


# Copy the flag versions
flag_dir = ms_active + ".flagversions"
flag_out_dir = os.path.join(out_path, os.path.basename(flag_dir))
os.mkdir(flag_out_dir)
copy_tree(flag_dir, flag_out_dir)

# Copy the cal tables
cal_dir = os.path.join(ms_path, "final_caltables")
cal_out_dir = os.path.join(out_path, os.path.basename(cal_dir))
os.mkdir(cal_out_dir)
copy_tree(cal_dir, cal_out_dir)

# Reduction logs
logs_dir = os.path.join(ms_path, "logs")
# Copy any extra log files into the logs dir.
for log_file in glob(ms_path + "/*.log"):
    shutil.move(log_file, logs_dir)
logs_out_dir = os.path.join(out_path, os.path.basename(logs_dir))
os.mkdir(logs_out_dir)
copy_tree(logs_dir, logs_out_dir)

# Test images
img_dir = os.path.join(ms_path, "test_images")
img_out_dir = os.path.join(out_path, os.path.basename(img_dir))
os.mkdir(img_out_dir)
copy_tree(img_dir, img_out_dir)

# Copy the reduction pipeline plots/outputs
weblog_file = os.path.join(ms_path, "weblog.speclines.tgz")
copy_file(weblog_file, out_path)

# Copy the pipeline.restore file
restore_file = os.path.join(ms_path, "pipeline_shelf.restore")
copy_file(restore_file, out_path)

# List the files in out_path
os.listdir(out_path)

# Now go and zip the output path
out_path_parent = os.path.dirname(out_path)
os.chdir(out_path_parent)
os.system("tar -zcf {0}.tar.gz {0}".format(os.path.basename(out_path)))
