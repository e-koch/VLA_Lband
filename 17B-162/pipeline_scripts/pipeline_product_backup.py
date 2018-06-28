
'''
Backup an MS's flag backup and the final calibration tables.

Works with VLA pipeline with CASA 5.1.2
'''

import sys
import os
from distutils.dir_util import copy_tree
from distutils.file_util import copy_file
from glob import glob

ms_path = sys.argv[1]
out_path = sys.argv[2]
cleanup = True if sys.argv[3] == "T" else False

print(ms_path, out_path)

# Continuum or lines?
if "speclines" in ms_path:
    track_type = 'speclines'
elif 'continuum' in ms_path:
    track_type = 'continuum'
else:
    raise ValueError("Cannot disinguish if this is a line or continuum ms.")

ms_path_join = lambda x: os.path.join(ms_path, x)
out_path_join = lambda x: os.path.join(out_path, x)

if not os.path.exists(out_path):
    os.mkdir(out_path)

ms_active = glob(ms_path_join("*.{}.ms".format(track_type)))

# There should only be one...
assert len(ms_active) == 1

ms_active = ms_active[0]

print(ms_active)

# Copy the cal tables and flag versions
files = glob("{}.*".format(ms_active))
for f in files:
    if os.path.isdir(f):
        out_dir = out_path_join(os.path.basename(f))
        os.mkdir(out_dir)
        copy_tree(f, out_dir)
    else:
        copy_file(f, out_path)

# Also the flux gain cal
flux_cal = ms_path_join("fluxgaincal.g")
out_dir = out_path_join("fluxgaincal.g")
os.mkdir(out_dir)
copy_tree(flux_cal, out_dir)

# Copy the pipeline outputs
files = glob(ms_path_join("pipeline*"))
for f in files:
    if os.path.isdir(f):
        out_dir = out_path_join(os.path.basename(f))
        os.mkdir(out_dir)
        copy_tree(f, out_dir)
    else:
        copy_file(f, out_path)

# List the files in out_path
os.listdir(out_path)

# Remove unneeded pipeline products
if cleanup:

    os.system("rm -r {}/BPcal.b".format(ms_path))
    os.system("rm -r {}/BPinitialgain.g".format(ms_path))
    os.system("rm -r {}/calibrators.ms".format(ms_path))
    os.system("rm -r {}/delay.k".format(ms_path))
    os.system("rm -r {}/finalcalibrators.ms".format(ms_path))
    os.system("rm -r {}/fluxgaincalFcal.g".format(ms_path))
    os.system("rm -r {}/fluxphaseshortgaincal.g".format(ms_path))
    os.system("rm -r {}/scratch.*".format(ms_path))
    os.system("rm -r {}/semiFinaldelayinitialgain.g".format(ms_path))
    os.system("rm -r {}/testBPcal.b".format(ms_path))
    os.system("rm -r {}/testBPdinitialgain.g".format(ms_path))
    os.system("rm -r {}/testdelayinitialgain.g".format(ms_path))
    os.system("rm -r {}/testdelay.k".format(ms_path))
    os.system("rm -r {}/testgaincal.g".format(ms_path))

    # os.system("rm -r {}/image_outputs".format(ms_path))
    # os.system("rm -r {}/scan_plots".format(ms_path))

# Now go and zip the output path
out_path_parent = os.path.dirname(out_path)
os.chdir(out_path_parent)
os.system("tar -zcf {0}.tar.gz {0}".format(os.path.basename(out_path)))
