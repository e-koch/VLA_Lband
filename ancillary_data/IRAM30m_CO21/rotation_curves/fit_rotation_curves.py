
'''
Create moment maps compatible with DiskFit, run DiskFit, and create useful
outputs.
'''

from astropy.io import fits
from astropy import log
import numpy as np
import os
from threading import Thread
from galaxies import Galaxy

from cube_analysis.rotation_curves import run_diskfit

from paths import (iram_co21_data_path, iram_co21_analysispath)


# While DiskFit 1.2.1 should handle 64-bit arrays, I'm getting errors using
# the default dtype (>f4) for this data set. I think this is somehow a
# python-to-FITS issue. Create a few 32-bit moments to use instead.

out_path = iram_co21_data_path("moments_for_diskfit", no_check=True)

if not os.path.exists(out_path):
    os.mkdir(out_path)

mom1 = fits.open(iram_co21_data_path("m33.co21_iram.mom1.fits"))[0]
mom1.header["BITPIX"] = -32
mom1.data = mom1.data.astype(np.float32)

mom1_name = os.path.split(iram_co21_data_path("m33.co21_iram.mom1.fits"))[-1]

mom1.writeto(iram_co21_data_path("moments_for_diskfit/{}".format(mom1_name),
                                 no_check=True),
             overwrite=True)

peakvel = fits.open(iram_co21_data_path("m33.co21_iram.peakvels.fits"))[0]
peakvel.header["BITPIX"] = -32
peakvel.data = peakvel.data.astype(np.float32)

peakvel_name = os.path.split(iram_co21_data_path("m33.co21_iram.peakvels.fits"))[-1]

peakvel.writeto(iram_co21_data_path("moments_for_diskfit/{}".format(peakvel_name),
                                    no_check=True),
                overwrite=True)

gal = Galaxy("M33")

param_file = iram_co21_analysispath("rotation_curves/diskfit_params_moment1.inp")
data_path = iram_co21_data_path("", no_check=True)
fits_for_wcs = iram_co21_data_path("m33.co21_iram.mom1.fits")

thr1 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs),
              kwargs={"fit_model": True, "gal": gal})
log.info("Starting Centroid DiskFit run")
thr1.start()


param_file = iram_co21_analysispath("rotation_curves/diskfit_params_peakvels.inp")
data_path = iram_co21_data_path("", no_check=True)
fits_for_wcs = iram_co21_data_path("m33.co21_iram.peakvels.fits")

thr2 = Thread(target=run_diskfit, args=(param_file, data_path, fits_for_wcs),
              kwargs={"fit_model": True, "gal": gal})
log.info("Starting Peak Velocity DiskFit run")
thr2.start()

thr1.join()
log.info("Finished Centroid DiskFit run")
thr2.join()
log.info("Finished Peak Velocity DiskFit run")
