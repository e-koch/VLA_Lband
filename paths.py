import os
import socket
from functools import partial


'''
Common set of paths giving the location of data products.
'''


def name_return_check(filename, path, no_check=False):
    full_path = os.path.join(path, filename)

    if not os.path.exists(full_path) and not no_check:
        raise OSError("{} does not exist.".format(full_path))

    return full_path


if socket.gethostname() == 'ewk':
    root = os.path.expanduser('~/Dropbox/code_development/VLA_Lband/')
    data_path = "/mnt/MyRAID/M33/"
# Add in path for NRAO and cloud instances
elif socket.gethostname() == 'caterwauler':
    root = os.path.expanduser('~/Dropbox/code_development/VLA_Lband/')
    data_path = os.path.expanduser("~/volume/data/")
# NRAO
elif "nmpost" in socket.gethostname():
    root = os.path.expanduser("~/VLA_Lband")
    data_path = os.path.expanduser("~/data")
elif "segfault" == socket.gethostname():
    root = os.path.expanduser("~/Dropbox/code_development/VLA_Lband/")
    data_path = "/mnt/bigdata/ekoch/M33"

c_path = os.path.join(root, '14B-088')
archival_path = os.path.join(root, 'AT0206')
a_path = os.path.join(root, '16B')
archival_12_path = os.path.join(root, '12A-403')

c_hi_analysispath = \
    partial(name_return_check, path=os.path.join(c_path, 'HI/analysis'))
archival_hi_analysispath = os.path.join(archival_path, 'AT0206/Analysis')

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
    os.path.join(data_path, "14B-088/HI/full_imaging_wGBT/")
# fourteenB_HI_data_path = os.path.join(data_path, "14B-088/HI/full_imaging/")
arecibo_HI_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "Arecibo/"))

iram_co21_data_path = partial(name_return_check,
                              path=os.path.join(data_path, "co21/"))

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
