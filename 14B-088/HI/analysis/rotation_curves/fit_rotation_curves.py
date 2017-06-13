
'''
Create moment maps compatible with DiskFit, run DiskFit, and create useful
outputs.
'''

from astropy.io import fits
import numpy as np
import os
from threading import Thread

from cube_analysis.rotation_curves import run_diskfit

from paths import (fourteenB_HI_file_dict, fourteenB_HI_data_path,
                   c_hi_analysispath)


# While DiskFit 1.2.1 should handle 64-bit arrays, I'm getting errors using
# the default dtype (>f4) for this data set. I think this is somehow a
# python-to-FITS issue. Create a few 32-bit moments to use instead.

out_path = fourteenB_HI_data_path("moments_for_diskfit", no_check=True)

if not os.path.exists(out_path):
    os.mkdir(out_path)

mom1 = fits.open(fourteenB_HI_file_dict["Moment1"])[0]
mom1.header["BITPIX"] = -32
mom1.data = mom1.data.astype(np.float32)

mom1_name = os.path.split(fourteenB_HI_file_dict["Moment1"])[-1]

mom1.write(fourteenB_HI_data_path("moments_for_diskfit/{}".format(mom1_name)),
           overwrite=True)

peakvel = fits.open(fourteenB_HI_file_dict["PeakVels"])[0]
peakvel.header["BITPIX"] = -32
peakvel.data = peakvel.data.astype(np.float32)

peakvel_name = os.path.split(fourteenB_HI_file_dict["PeakVels"])[-1]

peakvel.write(fourteenB_HI_data_path("moments_for_diskfit/{}".format(peakvel_name)),
              overwrite=True)

# Now run each model. Submit as a thread so we can run these simultaneously.

param_file = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_nowarp_noradial_noasymm.inp")
data_path = fourteenB_HI_data_path("", no_check=True)
fits_for_wcs = fourteenB_HI_file_dict["Moment1"]

thr1 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs))
thr1.start()


param_file = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_peakvels_nowarp_noradial_noasymm.inp")
data_path = fourteenB_HI_data_path("", no_check=True)
fits_for_wcs = fourteenB_HI_file_dict["PeakVels"]

thr2 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs))
thr2.start()


thr1.join()
thr2.join()