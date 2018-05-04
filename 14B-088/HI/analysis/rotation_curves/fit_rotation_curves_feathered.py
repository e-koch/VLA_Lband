
'''
Create moment maps compatible with DiskFit, run DiskFit, and create useful
outputs.
'''

from astropy.io import fits
import numpy as np
import os
from threading import Thread
from galaxies import Galaxy

from cube_analysis.rotation_curves import run_diskfit

from paths import (fourteenB_wGBT_HI_file_dict, fourteenB_HI_data_wGBT_path,
                   c_hi_analysispath)


# While DiskFit 1.2.1 should handle 64-bit arrays, I'm getting errors using
# the default dtype (>f4) for this data set. I think this is somehow a
# python-to-FITS issue. Create a few 32-bit moments to use instead.

out_path = fourteenB_HI_data_wGBT_path("moments_for_diskfit", no_check=True)

if not os.path.exists(out_path):
    os.mkdir(out_path)

mom1 = fits.open(fourteenB_wGBT_HI_file_dict["Moment1"])[0]
mom1.header["BITPIX"] = -32
mom1.data = mom1.data.astype(np.float32)

mom1_name = os.path.split(fourteenB_wGBT_HI_file_dict["Moment1"])[-1]

mom1.writeto(fourteenB_HI_data_wGBT_path("moments_for_diskfit/{}".format(mom1_name), no_check=True),
             overwrite=True)

peakvel = fits.open(fourteenB_wGBT_HI_file_dict["PeakVels"])[0]
peakvel.header["BITPIX"] = -32
peakvel.data = peakvel.data.astype(np.float32)

peakvel_name = os.path.split(fourteenB_wGBT_HI_file_dict["PeakVels"])[-1]

peakvel.writeto(fourteenB_HI_data_wGBT_path("moments_for_diskfit/{}".format(peakvel_name), no_check=True),
                overwrite=True)

# Now run each model. Submit as a thread so we can run these simultaneously.

gal = Galaxy("M33")

param_file = c_hi_analysispath("rotation_curves/diskfit_params_feathered/diskfit_params_nowarp_noradial_noasymm.inp")
data_path = fourteenB_HI_data_wGBT_path("", no_check=True)
fits_for_wcs = fourteenB_wGBT_HI_file_dict["Moment1"]

thr1 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs),
              kwargs={"fit_model": True, "gal": gal})
thr1.start()
thr1.join()

param_file = c_hi_analysispath("rotation_curves/diskfit_params_feathered/diskfit_params_peakvels_nowarp_noradial_noasymm.inp")
data_path = fourteenB_HI_data_wGBT_path("", no_check=True)
fits_for_wcs = fourteenB_wGBT_HI_file_dict["PeakVels"]

thr2 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs),
              kwargs={"fit_model": True, "gal": gal})
thr2.start()
thr2.join()

param_file = c_hi_analysispath("rotation_curves/diskfit_params_feathered/diskfit_params_peakvels_nowarp_radial_noasymm.inp")
data_path = fourteenB_HI_data_wGBT_path("", no_check=True)
fits_for_wcs = fourteenB_wGBT_HI_file_dict["PeakVels"]

thr3 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs),
              kwargs={"fit_model": True, "gal": gal})
thr3.start()
thr3.join()
