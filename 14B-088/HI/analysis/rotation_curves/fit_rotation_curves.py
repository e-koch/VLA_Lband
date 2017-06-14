
'''
Create moment maps compatible with DiskFit, run DiskFit, and create useful
outputs.
'''

from astropy.io import fits
from astropy import log
from spectral_cube import SpectralCube
import numpy as np
import os
from threading import Thread
from galaxies import Galaxy

from cube_analysis.rotation_curves import (run_diskfit, return_smooth_model,
                                           update_galaxy_params)
from cube_analysis.spectral_fitting import cube_fitter
from cube_analysis.spectral_fitting.gausshermite import herm_gauss_peak

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

mom1.writeto(fourteenB_HI_data_path("moments_for_diskfit/{}".format(mom1_name), no_check=True),
             overwrite=True)

peakvel = fits.open(fourteenB_HI_file_dict["PeakVels"])[0]
peakvel.header["BITPIX"] = -32
peakvel.data = peakvel.data.astype(np.float32)

peakvel_name = os.path.split(fourteenB_HI_file_dict["PeakVels"])[-1]

peakvel.writeto(fourteenB_HI_data_path("moments_for_diskfit/{}".format(peakvel_name), no_check=True),
                overwrite=True)

# For comparison, we want to fit to the peak from a Gaussian-Hermite fit.
# This is (currently) the only time this is used, so do the fitting here, save
# the peak position.
# log.info("Making Gaussian-Hermite Peak array")
# cube = SpectralCube.read(fourteenB_HI_file_dict["Cube"])
# herm_gauss_peak_array = cube_fitter(cube, herm_gauss_peak, verbose=True,
#                                     num_cores=1)[0].squeeze()

# ghpeak = fits.PrimaryHDU(herm_gauss_peak_array, header=mom1.header)
# ghpeak.header["BITPIX"] = -32
# ghpeak.data = ghpeak.data.astype(np.float32)

# ghpeak_name = ".".join((".".split(mom1_name)[:-2] + "gh_peak.fits"))

# ghpeak.writeto(fourteenB_HI_data_path("moments_for_diskfit/{}".format(ghpeak_name), no_check=True),
#                overwrite=True)

# Now run each model. Submit as a thread so we can run these simultaneously.

gal = Galaxy("M33")

param_file = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_nowarp_noradial_noasymm.inp")
data_path = fourteenB_HI_data_path("", no_check=True)
fits_for_wcs = fourteenB_HI_file_dict["Moment1"]

thr1 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs),
              kwargs={"fit_model": True, "gal": gal})
log.info("Starting Centroid DiskFit run")
thr1.start()


param_file = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_peakvels_nowarp_noradial_noasymm.inp")
data_path = fourteenB_HI_data_path("", no_check=True)
fits_for_wcs = fourteenB_HI_file_dict["PeakVels"]

thr2 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs),
              kwargs={"fit_model": True, "gal": gal})
log.info("Starting Peak Velocity DiskFit run")
thr2.start()

# param_file = c_hi_analysispath("rotation_curves/diskfit_params/diskfit_params_ghfit_nowarp_noradial_noasymm.inp")
# data_path = fourteenB_HI_data_path("", no_check=True)
# # This has the same WCS info.
# fits_for_wcs = fourteenB_HI_file_dict["PeakVels"]

# thr3 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs))
# log.info("Starting GH Peak DiskFit run")
# thr3.start()

thr1.join()
log.info("Finished Centroid DiskFit run")
thr2.join()
log.info("Finished Peak Velocity DiskFit run")
# thr3.join()
# log.info("Finished GH Peak DiskFit run")

# Now fit the rotation models to make a smooth model

