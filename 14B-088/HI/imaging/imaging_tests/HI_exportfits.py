
'''
Calculate common statistics in each of the test images, looking for any key
differences.
'''

import os
from glob import glob
import numpy as np

data_path = os.path.expanduser("~/MyRAID/M33/14B-088/HI/channel_testing")
# data_path = os.path.expanduser("/Volumes/Travel_Data/M33/testing")


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

for folder in folders:

    os.system("rm {}/*.fits".format(folder))

    # Now try opening the image, if it exists
    imagename = os.path.join(folder, os.path.basename(folder) + ".clean.image")
    feather_imagename = os.path.join(folder, os.path.basename(folder) +
                                     ".clean.image.feathered")
    if os.path.exists(feather_imagename):
        image = feather_imagename
        fitsimage = os.path.join(folder, feather_imagename + ".fits")

    elif os.path.exists(imagename):
        image = imagename
        fitsimage = os.path.join(folder, imagename + ".fits")
    else:
        # Cleaning failed. Skip.
        continue

    exportfits(imagename=image, fitsimage=fitsimage, velocity=True,
               history=False, dropdeg=True)
