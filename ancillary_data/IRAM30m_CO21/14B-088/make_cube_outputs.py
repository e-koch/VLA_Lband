
'''
Make a signal mask and moment maps for the CO(2-1) data cube.
'''

from cube_analysis import run_pipeline

from paths import iram_co21_data_path


# Find a signal mask and derive moment arrays
run_pipeline(iram_co21_data_path("m33.co21_iram.fits"),
             iram_co21_data_path("", no_check=True),
             masking_kwargs={"method": "ppv_dilation",
                             "save_cube": True,
                             "noise_map": iram_co21_data_path("m33.rms.fits"),
                             "min_sig": 3,
                             "max_sig": 5,
                             "min_pix": 27,
                             },
             moment_kwargs={"num_cores": 6,
                            "verbose": True})
