
'''
Given a imagename prefix and path, gather all saved channel outputs
and create a table.
'''

import sys
import os
from glob import glob
import numpy as np
from astropy.table import Table, Column


path = sys.argv[1]

# W/o _channel_{num} ending
filename = sys.argv[2]

# stage 1 or 2
stage = int(sys.argv[3])

if stage == 1:
    file_end = ""
elif stage == 2:
    file_end = "_stage2"
    file_end_stage1 = ""

else:
    raise ValueError("stage must be '1' or '2'.")

folders = glob(os.path.join(path, "channel_*"))

chan_nums = [int(folder.split("_")[-1]) for folder in folders]
chan_nums.sort()

for i in chan_nums:
    npy_file = glob(os.path.join(path, "channel_{}".format(i),
                                 "{0}_channel_{1}.results_dict{2}.npy"
                                 .format(filename, i, file_end)))
    if len(npy_file) == 1:
        out_dict = np.load(npy_file[0], allow_pickle=True, encoding='latin1').item()
        out_dict['stage'] = stage
    else:
        # If stage 2, check for a stage 1 file for those that didn't need
        # further cleaning
        if stage == 1:
            continue

        npy_file = glob(os.path.join(path, "channel_{}".format(i),
                                     "{0}_channel_{1}.results_dict{2}.npy"
                                     .format(filename, i, file_end_stage1)))
        if len(npy_file) == 0:
            continue

        # out_dict = np.load(npy_file[0]).item()
        # When loading in python3:
        out_dict = np.load(npy_file[0], allow_pickle=True, encoding='latin1').item()

        out_dict['stage'] = 1

    out_dict['channel'] = i

    # Setup a dictionary
    try:
        all_dict
    except NameError:
        all_dict = dict.fromkeys(out_dict)

        for key in out_dict:
            all_dict[key] = []

    for key in out_dict:
        all_dict[key].append(out_dict[key])

# for i, fil in enumerate(files):

#     out_dict = np.load(fil).item()

#     # Setup a dictionary
#     if i == 0:
#         all_dict = dict.fromkeys(out_dict)

#         for key in all_dict:
#             all_dict[key] = []

#     for key in out_dict:
#         all_dict[key].append(out_dict[key])


# Now convert the dictionary into an astropy table

table = Table()

for key in all_dict:

    # The shape will change with number of minor cycles
    # Need to handle this some other way...
    if key == 'summaryminor' or key == 'summarymajor':
        continue

    table[key] = Column(all_dict[key])


stop_codes = {0: "Not reached", 1: "Reached niter",
              2: "Reached threshold",
              3: "Stop flag", 4: "No change after cycle",
              5: "Diverging. Peak residual increase.",
              6: "Diverging. Peak min. residual increase.",
              7: "Empty clean mask",
              8: "Reached nsigma threshold"}

# Empty string array set to longest message in stop_codes
str_outs = np.empty(table['nsigma'].shape, dtype='S39')

str_outs = []

for code in table['stopcode']:

    if code < 0 or code > 8:
        raise ValueError("Found code {0} not within 0 and 8. Check"
                         " output.".format(code))

    str_outs.append(stop_codes[code])

table['stopcode_exp'] = Column(str_outs)

# Save the table in the given path

table.write(os.path.join(path, filename + file_end + ".h5"), path='data',
            format='hdf5', overwrite=True)
