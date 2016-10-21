
Steps in per-channel-imaging pipeline:
--------------------------------------

1. Transform MS into desired velocity frame and channel width. (e.g., run `mstransform`).
2. Split individual channels out of MS (can be time-consuming for large MS; `ms_channel_split.py`)
3. Regrid mask and model (if using) to the dimensions of the output cube.
4. Split out individual mask and model channels. (`cube_channel_split_template.py`)
5. Generate jobs for imaging each channel. (`job_generator_template.py`, `single_channel_clean_template.py`)
6. Search each CASA log to extract whether clean converged (plus number of iterations, time that it took to run). (`get_clean_results_template.py`)
7. Combine channel products into final cubes. (`create_cube_template.py`)
