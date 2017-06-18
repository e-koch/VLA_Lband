
'''
Make signal masks and compute the moments.
'''

from astropy import log
from cube_analysis import run_pipeline

from paths import (fourteenB_HI_file_dict, fourteenB_HI_data_path,
                   fourteenB_wGBT_HI_file_dict, fourteenB_HI_data_wGBT_path)


# VLA-only cube
log.info("Masking and moments for the VLA-only cube")
run_pipeline(fourteenB_HI_file_dict["Cube"],
             fourteenB_HI_data_path("", no_check=True),
             masking_kwargs={"method": "ppv_connectivity",
                             "save_cube": True,
                             "is_huge": True,
                             "noise_map": None,
                             "smooth_chans": 31,
                             "min_chan": 10,
                             "peak_snr": 5.,
                             "min_snr": 2,
                             "edge_thresh": 1,
                             },
             moment_kwargs={"num_cores": 6,
                            "verbose": True})

# VLA+GBT cube
log.info("Masking and moments for the VLA+GBT cube")
run_pipeline(fourteenB_wGBT_HI_file_dict["Cube"],
             fourteenB_HI_data_wGBT_path("", no_check=True),
             masking_kwargs={"method": "ppv_connectivity",
                             "save_cube": True,
                             "is_huge": True,
                             "noise_map": None,
                             "smooth_chans": 31,
                             "min_chan": 10,
                             "peak_snr": 5.,
                             "min_snr": 2,
                             "edge_thresh": 1,
                             },
             moment_kwargs={"num_cores": 6,
                            "verbose": True})
