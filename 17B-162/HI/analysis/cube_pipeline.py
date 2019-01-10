
'''
Make signal masks and compute the moments.
'''

from astropy import log
from cube_analysis import run_pipeline


run_02kms = True
run_1kms = True

num_cores = 4

if run_02kms:
    log.info("Running 0.21 km/s cubes")

    from paths import (seventeenB_HI_data_02kms_path,
                       seventeenB_HI_data_02kms_wGBT_path)

    # VLA-only cube
    log.info("Masking and moments for the VLA-only cube")
    run_pipeline(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.fits"),
                 seventeenB_HI_data_02kms_path("", no_check=True),
                 pb_file=seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.pb.fits"),
                 pb_lim=0.5, convolve_to_common_beam=True,
                 masking_kwargs={"method": "ppv_connectivity",
                                 "save_cube": True,
                                 "is_huge": True,
                                 "smooth_chans": 31,
                                 "min_chan": 10,
                                 "peak_snr": 4.,
                                 "min_snr": 2,
                                 "edge_thresh": 1,
                                 "pb_map_name": seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.pb.pbcov_gt_0.5_masked.fits")
                                 },
                 moment_kwargs={"num_cores": num_cores,
                                "verbose": True,
                                "chunk_size": 1e5,
                                "make_peakvels": False},
                 combeam_kwargs={'epsilon': 9e-4})

    # VLA+GBT cube
    log.info("Masking and moments for the VLA+GBT cube")
    run_pipeline(seventeenB_HI_data_02kms_wGBT_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.GBT_feathered.fits"),
                 seventeenB_HI_data_02kms_wGBT_path("", no_check=True),
                 pb_file=seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.pb.fits"),
                 pb_lim=0.5, convolve_to_common_beam=True,
                 masking_kwargs={"method": "ppv_connectivity",
                                 "save_cube": True,
                                 "is_huge": True,
                                 "smooth_chans": 31,
                                 "min_chan": 10,
                                 "peak_snr": 4.,
                                 "min_snr": 2,
                                 "edge_thresh": 1,
                                 "pb_map_name": seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.pb.pbcov_gt_0.5_masked.fits")
                                 },
                 moment_kwargs={"num_cores": num_cores,
                                "verbose": True,
                                "chunk_size": 1e5,
                                "make_peakvels": False},
                 combeam_kwargs={'epsilon': 9e-4})

if run_1kms:
    log.info("Running 1.0 km/s cubes")

    from paths import (seventeenB_HI_data_1kms_path,
                       seventeenB_HI_data_1kms_wGBT_path)

    # VLA-only cube
    log.info("Masking and moments for the VLA-only cube")
    run_pipeline(seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.image.pbcor.fits"),
                 seventeenB_HI_data_1kms_path("", no_check=True),
                 pb_file=seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.pb.fits"),
                 pb_lim=0.5, convolve_to_common_beam=True,
                 masking_kwargs={"method": "ppv_connectivity",
                                 "save_cube": True,
                                 "is_huge": True,
                                 "smooth_chans": 6,
                                 "min_chan": 4,
                                 "peak_snr": 4.,
                                 "min_snr": 2,
                                 "edge_thresh": 1,
                                 "verbose": True,
                                 },
                 moment_kwargs={"num_cores": num_cores,
                                "verbose": True,
                                "chunk_size": 1e5,
                                "make_peakvels": False},
                 combeam_kwargs={'epsilon': 9e-4})

    # VLA+GBT cube
    log.info("Masking and moments for the VLA+GBT cube")
    run_pipeline(seventeenB_HI_data_1kms_wGBT_path("M33_14B_17B_HI_contsub_width_1kms.image.pbcor.GBT_feathered.fits"),
                 seventeenB_HI_data_1kms_wGBT_path("", no_check=True),
                 pb_file=seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.pb.fits"),
                 pb_lim=0.5, convolve_to_common_beam=True,
                 masking_kwargs={"method": "ppv_connectivity",
                                 "save_cube": True,
                                 "is_huge": True,
                                 "smooth_chans": 6,
                                 "min_chan": 4,
                                 "peak_snr": 4.,
                                 "min_snr": 2,
                                 "edge_thresh": 1,
                                 "verbose": True,
                                 },
                 moment_kwargs={"num_cores": num_cores,
                                "verbose": True,
                                "chunk_size": 1e5,
                                "make_peakvels": False},
                 combeam_kwargs={'epsilon': 9e-4})
