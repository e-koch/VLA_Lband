Contains imaging scripts for line and continuum SPWs to be used in the VLA pipeline environment.

To use, add a symbolic link of this folder into the folder that contains the pipeline scripts.

To run as part of the pipeline, add an additional prompted argument in `EVLA_pipe_startup.py` to enable imaging at the end of the pipeline (`test_imaging` as a boolean). If enabled, also prompt at that time for the sources to image (e.g., `imaging_sources='M33, 3C48'`, where wildcards are added for imaging mosaics). Then, add a new line to `EVLA_pipeline.py` that runs one of these scripts if `test_imaging == True`. To restore the environment, add `test_imaging` and `imaging_souces` to the pipeline restore list.