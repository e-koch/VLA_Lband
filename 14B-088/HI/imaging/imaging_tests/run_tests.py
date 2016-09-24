
'''
Generate jobs on Jasper to run HI channel imaging tests.
'''

import os
from shutil import copytree
from copy import copy
import subprocess

# Job parameters
NODE = "1"
PROCS = "12"
PMEM = "4000mb"
HOURS = "72"

# Run imaging tests w/ different parameters/CASA versions
output_path = os.path.expanduser("~/m33/14B-088/testing")

os.chdir(output_path)

filename = os.path.expanduser("~/code_repos/VLA_Lband/14B-088/HI/imaging/"
                              "imaging_tests/HI_testing_channel_clean.py")

ms_name = "14B-088_HI_LSRK.ms.contsub_channel_1000.ms"
modelname = "M33_14B-088_HI_model_channel_330.image"
maskname = "M33_14B-088_HI_mask_modified_channel_330.image"

# Note that "fixspwbackport" must be run on the channel ms to work in 4.2.2
# https://safe.nrao.edu/wiki/bin/view/Software/CASAUserTestingMinutes20141021
# I did this by hand, not in the script

versions = [os.path.expanduser("~/casapy-42.2.30986-pipe-1-64b/bin/casa-4.2.2"),
            os.path.expanduser("~/casa-release-4.3.1-el6/bin/casa-4.3.1"),
            os.path.expanduser("~/casa-release-4.4.0-el6/bin/casa-4.4"),
            os.path.expanduser("~/casa-release-4.5.3-el6/bin/casa-4.5.3"),
            os.path.expanduser("~/casa-release-4.6.0-el6/bin/casa")]

call = 'qsub -N JOB_NAME -l nodes=NODE:ppn=PROCS,pmem=PMEM,' \
    'walltime=HOURS:00:00 -d . <<< "VERSION --logfile JOB_NAME.log -c' \
    ' FILENAME MS_NAME MODEL MASK FIELD MSCALE TCLEAN"'

# Set job parameters. Name is done in the loop.
call = call.replace("NODE", NODE)
call = call.replace("PROCS", PROCS)
call = call.replace("PMEM", PMEM)
call = call.replace("HOURS", HOURS)

call = call.replace("FILENAME", filename)
call = call.replace("MS_NAME", ms_name)

# Now loop through all combinations
for version in versions:
    for tclean in ["T", "F"]:

        if os.path.basename(version) == "casa":
            major, minor, revision = 4, 6, 0
        else:
            try:
                major, minor, revision = \
                    os.path.basename(version).split("-")[1].split('.')
            except ValueError:
                major, minor = \
                    os.path.basename(version).split("-")[1].split('.')
                revision = 0

        casa_version = 100 * int(major) + 10 * int(minor) + int(revision)

        # Skip if on a version w/o tclean
        if tclean == "T" and int(minor) < 5:
            continue

        for model in [None, modelname]:
            for mask in [None, maskname]:
                for fields in ["T", "F"]:
                    for mscale in ["T", "F"]:
                        JOB_NAME = "{0}.CASAVer_{1}.Model_{2}.Mask_{3}." \
                                   "AllFields_{4}.MScale_{5}" \
                                   ".Tclean_{6}".format(ms_name[:-3],
                                                        casa_version,
                                                        "T" if model is
                                                        not None else "F",
                                                        "T" if mask is
                                                        not None else "F",
                                                        fields,
                                                        mscale,
                                                        tclean)

                        new_call = copy(call)
                        new_call = new_call.replace("JOB_NAME", JOB_NAME)
                        new_call = new_call.replace("VERSION", version)
                        new_call = new_call.replace("TCLEAN", tclean)
                        new_call = new_call.replace("MODEL", str(model))
                        new_call = new_call.replace("MASK", str(mask))
                        new_call = new_call.replace("FIELD", fields)
                        new_call = new_call.replace("MSCALE", mscale)

                        job_folder = os.path.join(output_path, JOB_NAME)
                        # Assume that the data are already copied if the folder
                        # exists.
                        if not os.path.exists(job_folder):
                            os.mkdir(job_folder)

                            # Copy the ms, and the model and mask, if needed.
                            copytree(ms_name, os.path.join(job_folder,
                                                           ms_name))
                            if model is not None:
                                copytree(model, os.path.join(job_folder,
                                                             model))
                            if mask is not None:
                                copytree(mask, os.path.join(job_folder, mask))

                        os.chdir(job_folder)

                        sp = subprocess.Popen(["/bin/bash", "-i", "-c",
                                               new_call])
                        sp.communicate()

                        os.chdir(output_path)
