import os
import socket

if socket.gethostname() == 'ewk':
    root = os.path.expanduser('~/Dropbox/code_development/VLA_Lband/')
    data_path = "/media/eric/MyRAID/"
# Add in path for NRAO and cloud instances
elif socket.gethostname() == 'caterwauler':
    root = os.path.expanduser('~/Dropbox/code_development/VLA_Lband/')
    data_path = os.path.expanduser("~/volume/data/")
# NRAO
elif "nmpost" in socket.gethostname():
    root = os.path.expanduser("~/")
    data_path = os.path.expanduser("~/data")
# else:
#     root = os.path.expanduser('~/work/w51/alma/')

c_path = os.path.join(root, '14B-088')
archival_path = os.path.join(root, 'AT0206')
a_path = os.path.join(root, '16B')
archival_12_path = os.path.join(root, '12A-403')

c_hi_analysispath = os.path.join(c_path, 'HI/analysis')
archival_hi_analysispath = os.path.join(archival_path, 'AT0206/Analysis')

# Pipeline paths
fourteenB_pipe_path = os.path.join(c_path, "Cal_Scripts/EVLA_pipeline1.3.0") + "/"
sixteenB_pipe_path = os.path.join(a_path, "pipeline4.6.0_custom") + "/"
twelveA_pipe_path = os.path.join(archival_12_path, "pipeline4.6.0") + "/"

# Paths to common modules
image_script_path = os.path.join(root, 'imaging_pipeline')


def path(x, basepath):
    return os.path.join(basepath, x)
