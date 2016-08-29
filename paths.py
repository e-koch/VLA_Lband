import os
import socket

if socket.gethostname() == 'ewk':
    root = os.path.expanduser('~/Dropbox/code_development/VLA_Lband/')
    data_path = "/media/eric/MyRAID/"
# Add in path for NRAO and cloud instances
# elif socket.gethostname() == 'orion':
#     root = '/diskb//aginsbur/w51/merge/'
# else:
#     root = os.path.expanduser('~/work/w51/alma/')

c_path = os.path.join(root, '14B-088')
archival_path = os.path.join(root, 'AT0206')
a_path = os.path.join(root, '16B')

c_hi_analysispath = os.path.join(c_path, 'HI/analysis')
archival_hi_analysispath = os.path.join(archival_path, 'AT0206/Analysis')

# Pipeline paths
fourteenB_pipe_path = os.path.join(c_path, "Cal_Scripts/EVLA_pipeline1.3.0")
sixteenB_pipe_path = os.path.join(a_path, "pipeline4.6.0_custom")


def path(x, basepath):
    return os.path.join(basepath, x)
