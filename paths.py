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
# else:
#     root = os.path.expanduser('~/work/w51/alma/')

c_path = os.path.join(root, '14B-088')
archival_path = os.path.join(root, 'AT0206')
a_path = os.path.join(root, '16B')
archival_12_path = os.path.join(root, '12A-403')

c_hi_analysispath = \
    partial(name_return_check, path=os.path.join(c_path, 'HI/analysis'))
archival_hi_analysispath = os.path.join(archival_path, 'AT0206/Analysis')

# Pipeline paths
fourteenB_pipe_path = os.path.join(c_path, "Cal_Scripts/EVLA_pipeline1.3.0") + "/"
sixteenB_pipe_path = os.path.join(a_path, "pipeline4.6.0_custom") + "/"
twelveA_pipe_path = os.path.join(archival_12_path, "pipeline4.6.0") + "/"

# Paths to common modules
image_script_path = os.path.join(root, 'imaging_pipeline')

# Data paths
fourteenB_HI_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "VLA/14B-088/HI/full_imaging_1016/"))
# fourteenB_HI_data_path = os.path.join(data_path, "14B-088/HI/full_imaging/")
arecibo_HI_data_path = \
    partial(name_return_check,
            path=os.path.join(data_path, "Arecibo/"))

# Paper figures path
papers_path = os.path.expanduser("~/Dropbox/My_Papers/")
paper1_figures_path = \
    os.path.join(papers_path, "In Prep/m33_14b088_hi/figures/")
