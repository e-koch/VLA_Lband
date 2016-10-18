Custom Pipeline Changes
=======================

* `EVLA_pipeline.py` - only runs startup and deterministic flagging. The MS is then split into continuum and line SPWs, and custom pipelines for each are then completed.
* `EVLA_pipe_startup.py` - extra prompts for enabling test imaging.
* `EVLA_pipe_mixed_setup_split.py` - splits an MS into continuum and line spectral windows.
* `EVLA_pipeline_continuum.py` - runs the standard pipeline on the continuum MS
* `EVLA_pipeline_lines.py` - runs a modified pipeline on the line MS. `rflag` is not run on the target fields and `statwt` uses only a portion of the channels to calculate weights, making the assumption that the line is confined to the center third of the SPW.
* `EVLA_pipeline_statwt_lines.py` & `EVLA_pipeline_targetflag_lines.py` - line specific scripts. See directly above.
* `EVLA_pipeline_statwt_lines.py` - Calculates weights based on assumption that line emission is confined to the 40% of the band. To avoid the band edges, weights are calculated from the 10-30th and 70-90th channel percentiles in the SPW.Pqui
* The `pipeline_path` is set by `paths.py` at the top of the repository (and locally linked here). Before running the pipeline, add `sys.path.append("PATH_TO_VLA_Lband")` in CASA.
* `EVLA_pipe_testimage_cont.py` - Create dirty images of the calibrator, target, or both.
* `EVLA_pipe_testimage_lines.py` - Same as above. Images a single wide channel (5x spectral resolution) at the middle of the band.
* `EVLA_pipe_restore.list` - add parameters for test imaging
* `EVLA_pipe_restart_lines.py` - Changed some script names to use the line versions.
* Both restart scripts - Added test imaging scripts to the end. These are skipped if the `test_imaging` is `False`.
* Added `EVLA_pipe_fake_flagall.py` when running the split MS pipelines. This calculates the input flag statistics only, since all deterministic flagging is performed prior to splitting. Setting these parameters is necessary for other portions of the pipeline, though.

Project Specific Changes
========================

* Fields for polarization calibration are dropped from the line MS (`EVLA_pipe_mixed_setup_split.py`)
* 
