
'''
Compare feathering with different SD data.
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
                              "imaging_tests/sd_combo_tests/"
                              "HI_single_channel_tclean.py")

ms_name = "14B-088_HI_LSRK.ms.contsub_channel_1000.ms"
maskname = "M33_14B-088_HI_mask_modified_channel_330.image"

modelnames = {"Arecibo": "M33_14B-088_HI_model_channel_330.image",
              "Effelsburg": "",
              "GBT": "",
              "None": "None"}


casa_call = os.path.expanduser("~/casa-release-4.7.0-el6/bin/casa")

call = 'qsub -N JOB_NAME -l nodes=NODE:ppn=PROCS,pmem=PMEM,' \
    'walltime=HOURS:00:00 -d . <<< "VERSION --logfile JOB_NAME.log -c' \
    'FILENAME MS_NAME MODEL USE_MODEL MASK OUT_ROOT"'

# Set job parameters. Name is done in the loop.
call = call.replace("NODE", NODE)
call = call.replace("PROCS", PROCS)
call = call.replace("PMEM", PMEM)
call = call.replace("HOURS", HOURS)

call = call.replace("FILENAME", filename)
call = call.replace("MS_NAME", ms_name)

# Now loop through all combinations

for tele in modelnames:
    for use_cleanmodel in [True, False]:
        # Don't run the no-SD case twice.
        if tele == "None" and not use_cleanmodel:
            continue

        JOB_NAME = "{0}.Model_{1}.AsCleanModel_{2}"\
            .format(ms_name[:-3],
                    tele,
                    use_cleanmodel)

        new_call = copy(call)
        new_call = new_call.replace("JOB_NAME", JOB_NAME)
        new_call = new_call.replace("MODEL", modelnames[tele])
        new_call = new_call.replace("USE_MODEL", use_cleanmodel)
        new_call = new_call.replace("MASK", maskname)
        new_call = new_call.replace("OUT_ROOT", JOB_NAME)

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
