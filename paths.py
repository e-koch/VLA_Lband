import os
import socket
from functools import partial
import glob

'''
Common set of paths giving the location of data products.
'''


def name_return_check(filename, path, no_check=False):
    full_path = os.path.join(path, filename)

    if not os.path.exists(full_path) and not no_check:
        raise OSError("{} does not exist.".format(full_path))

    return full_path


if socket.gethostname() == 'ewk':
    root = os.path.expanduser('~/ownCloud/code_development/VLA_Lband/')
    # data_path = "/mnt/MyRAID/M33/"
    data_path = "/home/eric/bigdata/ekoch/M33/"
# Add in path for NRAO and cloud instances
elif socket.gethostname() == 'caterwauler':
    root = os.path.expanduser('~/ownCloud/code_development/VLA_Lband/')
    data_path = os.path.expanduser("~/volume/data/")
# NRAO
elif "nmpost" in socket.gethostname():
    root = os.path.expanduser("~/VLA_Lband")
    data_path = os.path.expanduser("~/data")
elif "segfault" == socket.gethostname():
    root = os.path.expanduser("~/ownCloud/code_development/VLA_Lband/")
    data_path = "/mnt/bigdata/ekoch/M33"
elif "cedar.computecanada" in socket.gethostname():
    root = "/home/ekoch/code/VLA_Lband/"
    data_path = "/home/ekoch/project/ekoch/"
elif 'ewk-laptop' in socket.gethostname():
    root = os.path.expanduser("~/ownCloud/code_development/VLA_Lband/")
    data_path = os.path.expanduser("~/storage/M33")

c_path = os.path.join(root, '14B-088')
archival_path = os.path.join(root, 'AT0206')
a_path = os.path.join(root, '16B')
archival_12_path = os.path.join(root, '12A-403')
ancillary_path = os.path.join(root, 'ancillary_data')
b_path = os.path.join(root, '17B-162')

c_hi_analysispath = \
    partial(name_return_check, path=os.path.join(c_path, 'HI/analysis'))
archival_hi_analysispath = os.path.join(archival_path, 'AT0206/Analysis')

iram_co21_analysispath = \
    partial(name_return_check,
            path=os.path.join(ancillary_path, 'IRAM30m_CO21'))

# Pipeline paths
fourteenB_pipe_path = os.path.join(c_path, "Cal_Scripts/EVLA_pipeline1.3.0") + "/"
sixteenB_pipe_path = os.path.join(a_path, "pipeline4.7.1_custom") + "/"
twelveA_pipe_path = os.path.join(archival_12_path, "pipeline4.6.0") + "/"

# Paths to common modules
image_script_path = os.path.join(root, 'imaging_pipeline')

# Data paths
fourteenB_HI_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "VLA/14B-088/HI/full_imaging_noSD/"))
fourteenB_HI_data_wGBT_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "VLA/14B-088/HI/full_imaging_wGBT/"))

seventeenB_HI_data_02kms_path = \
    partial(name_return_check,
            path=os.path.join(data_path,
                              "VLA/17B-162/HI/full_imaging_02kms_noSD/"))
seventeenB_HI_data_02kms_wGBT_path = \
    partial(name_return_check,
            path=os.path.join(data_path,
                              "VLA/17B-162/HI/full_imaging_02kms_wGBT/"))

seventeenB_HI_data_1kms_path = \
    partial(name_return_check,
            path=os.path.join(data_path,
                              "VLA/17B-162/HI/full_imaging_1kms_noSD/"))
seventeenB_HI_data_1kms_wGBT_path = \
    partial(name_return_check,
            path=os.path.join(data_path,
                              "VLA/17B-162/HI/full_imaging_1kms_wGBT/"))

arecibo_HI_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "Arecibo/"))

ebhis_HI_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "EBHIS/"))

gbt_HI_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "GBT/"))

iram_co21_data_path = partial(name_return_check,
                              path=os.path.join(data_path, "co21/"))
iram_co21_14B088_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "co21/14B-088/"))

hst_phat_path = \
  partial(name_return_check,
          path=os.path.join(data_path, "Hubble/big/"))

# Paper figures path
papers_path = os.path.expanduser("~/Dropbox/My_Papers/")
paper1_figures_path = \
    lambda x: os.path.join(papers_path, "In Prep/m33-HI-paper1/figures/", x)
paper1_tables_path = \
    lambda x: os.path.join(papers_path, "In Prep/m33-HI-paper1/tables/", x)

# Proposal Figures
varfig_path = os.path.expanduser("~/Dropbox/Various Plots/Proposals")
proposal_figures_path = lambda x: os.path.join(varfig_path, x)

# All figures
fig_path = os.path.expanduser("~/Dropbox/Various Plots/M33/")
allfigs_path = lambda x: os.path.join(fig_path, x)
alltables_path = lambda x: os.path.join(fig_path, "tables", x)


def find_dataproduct_names(path):
    '''
    Given a path, return a dictionary of the data products with the name
    convention used in this repository.
    '''

    search_dict = {"Moment0": "mom0",
                   "Moment1": "mom1",
                   "LWidth": "lwidth",
                   "Skewness": "skewness",
                   "Kurtosis": "kurtosis",
                   "PeakTemp": "peaktemps",
                   "PeakVels": "peakvels.",
                   "Cube": "masked.fits",
                   "Source_Mask": "masked_source_mask.fits",
                   "CentSub_Cube": "masked.centroid_corrected",
                   "CentSub_Mask": "masked_source_mask.centroid_corrected",
                   "RotSub_Cube": "masked.rotation_corrected",
                   "RotSub_Mask": "masked_source_mask.rotation_corrected",
                   "PeakSub_Cube": "masked.peakvels_corrected",
                   "PeakSub_Mask": "masked_source_mask.peakvels_corrected"}

    found_dict = {}

    for filename in glob.glob(os.path.join(path, "*.fits")):

        for key in search_dict:
            if search_dict[key] in filename:
                found_dict[key] = filename
                search_dict.pop(key)
                break

    return found_dict


# Return dictionaries with names for the existing directories
fourteenB_HI_file_dict = \
    find_dataproduct_names(fourteenB_HI_data_path("", no_check=True))
fourteenB_wGBT_HI_file_dict = \
    find_dataproduct_names(fourteenB_HI_data_wGBT_path("", no_check=True))

seventeenB_02kms_HI_file_dict = \
    find_dataproduct_names(seventeenB_HI_data_02kms_path("", no_check=True))
seventeenB_02kms_wGBT_HI_file_dict = \
    find_dataproduct_names(seventeenB_HI_data_02kms_wGBT_path("",
                                                              no_check=True))

seventeenB_1kms_HI_file_dict = \
    find_dataproduct_names(seventeenB_HI_data_1kms_path("", no_check=True))
seventeenB_1kms_wGBT_HI_file_dict = \
    find_dataproduct_names(seventeenB_HI_data_1kms_wGBT_path("",
                                                             no_check=True))

if __name__ == "__main__":

    # Append the repo directory to the path so paths is importable
    os.sys.path.append(root)
