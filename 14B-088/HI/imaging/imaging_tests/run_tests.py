
'''
Generate jobs on Jasper to run HI channel imaging tests.
'''

import os
from shutil import copytree

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

versions = ["casa-4.2.2", "casa-4.3.1", "casa-4.4", "casa-4.5.3", "casa-4.6"]

call = 'qsub -N JOB_NAME -l nodes=NODE:ppn=PROCS,pmem=PMEM,' \
    'walltime=HOUR:00:00 -d . <<< "VERSION --logfile JOB_NAME.log -c' \
    ' FILENAME MS_NAME MODEL MASK FIELD MSCALE TCLEAN"'

# Set job parameters. Name is done in the loop.
call = call.replace("NODE", NODE)
call = call.replace("PROCS", PROCS)
call = call.replace("PMEM", PMEM)
call = call.replace("HOURS", HOURS)

# Now loop through all combinations
for version in versions:
    for tclean in ["T", "F"]:
        major, minor, revision = casadef.casa_version.split('.')
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

                        new_call = call.copy()
                        new_call = new_call.replace("JOB_NAME", JOB_NAME)
                        new_call = new_call.replace("VERSION", version)
                        new_call = new_call.replace("TCLEAN", tclean)
                        new_call = new_call.replace("MODEL", model)
                        new_call = new_call.replace("MASK", mask)
                        new_call = new_call.replace("FIELD", fields)
                        new_call = new_call.replace("MSCALE", mscale)

                        job_folder = os.path.join(output_path, JOB_NAME)
                        if not os.path.exists(job_folder):
                            os.mkdir(job_folder)

                        # Copy the ms, and the model and mask, if needed.
                        copytree(ms_name, os.path.join(job_folder, ms_name))
                        if model is not None:
                            copytree(model, os.path.join(job_folder, model))
                        if mask is not None:
                            copytree(mask, os.path.join(job_folder, mask))

                        os.chdir(job_folder)

                        os.system("'source .bashrc; shopt -s"
                                  " expand_aliases;" + call)

                        os.chdir(output_path)
