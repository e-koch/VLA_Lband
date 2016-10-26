
'''
Calculate common statistics in each of the test images, looking for any key
differences.
'''

import os
from glob import glob
import numpy as np

from astropy.table import Table
from astropy.utils.console import ProgressBar

data_path = os.path.expanduser("~/MyRAID/M33/14B-088/HI/imaging/testing")

def append_path(path):
    return os.path.join(data_path, path)

parameters = ["CASAVer", "Model", "Mask", "AllFields", "MScale", "Tclean"]

# Load in the mask
ia.open(append_path("M33_14B-088_HI_mask_channel_330.image"))
mask = ia.torecord()["imagearray"].squeeze() > 0
ia.close()
# For some reason, both axes need to be reversed
mask = mask[::-1, ::-1]


folders = glob(append_path("14B-088_HI_LSRK.ms.contsub_channel_1000.CASAVer*"))

parameter_array = np.zeros((len(parameters), len(folders)))

parameter_table = dict.fromkeys(parameters)

for key in parameter_table:
    parameter_table[key] = []

# Values to extract out.
parameter_table["std"] = []
parameter_table["sum"] = []
parameter_table["median"] = []
parameter_table["peak_res"] = []

for folder in ProgressBar(folders):

    # Slice out the first 3 since this is just the name
    split_name = folder.split(".")[3:]

    for name in split_name:
        for param in parameters[1:]:
            if "CASAVer" in name:
                parameter_table["CASAVer"].append(int(name.split("_")[-1]))
                break
            if param in name:
                parameter_table[param].append(1 if name.split("_")[-1] == "T" else 0)
                break

    # Now try opening the image, if it exists
    imagename = os.path.join(folder, os.path.basename(folder) + ".clean.image")
    if not os.path.exists(imagename):
        # print("No image for {}".format(folder))
        parameter_table["std"].append(np.NaN)
        parameter_table["sum"].append(np.NaN)
        parameter_table["median"].append(np.NaN)
        parameter_table["peak_res"].append(np.NaN)
        continue

    ia.open(imagename)
    arr = ia.torecord()["imagearray"].squeeze()
    arr[~mask] = np.NaN
    ia.close()

    # The numpy version in CASA doesn't have the nan version?
    parameter_table["std"].append(np.std(arr[np.isfinite(arr)]))
    parameter_table["sum"].append(np.sum(arr[np.isfinite(arr)]))
    parameter_table["median"].append(np.median(arr[np.isfinite(arr)]))

    # Also pick out the peak residual in the residual image.
    imagename = os.path.join(folder, os.path.basename(folder) + ".clean.residual")
    ia.open(imagename)
    arr = ia.torecord()["imagearray"].squeeze()
    arr[~mask] = np.NaN
    ia.close()

    parameter_table["peak_res"].append(np.max(arr[np.isfinite(arr)]))

t = Table(parameter_table)
t.write(append_path("property_values.csv"), format='ascii.csv')
