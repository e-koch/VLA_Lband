Custom Pipeline Changes
=======================

* `EVLA_pipeline.py` - only runs startup and deterministic flagging. The MS is then split into continuum and line SPWs, and custom pipelines for each are then completed.
* `EVLA_pipe_startup.py` - extra prompts for enabling test imaging.
* `EVLA_pipe_mixed_setup_split.py` - splits an MS into continuum and line spectral windows.
* `EVLA_pipeline_continuum.py` - runs the standard pipeline on the continuum MS
* `EVLA_pipeline_lines.py` - runs a modified pipeline on the line MS. `rflag` is not run on the target fields and `statwt` uses only a portion of the channels to calculate weights, making the assumption that the line is confined to the center third of the SPW.
* `EVLA_pipeline_statwt_lines.py` & `EVLA_pipeline_targetflag_lines.py` - line specific scripts. See directly above.
* The `pipeline_path` is set by `paths.py` at the top of the repository (and locally linked here). Before running the pipeline, add `sys.path.append("PATH_TO_VLA_Lband")` in CASA.
* `EVLA_pipe_testimage_cont.py` - Create dirty images of the calibrator, target, or both.
* `EVLA_pipe_testimage_lines.py` - Same as above. Images a single wide channel (5x spectral resolution) at the middle of the band.

Project Specific Changes
========================

* Fields for polarization calibration are dropped from the line MS (`EVLA_pipe_mixed_setup_split.py`)
* 
