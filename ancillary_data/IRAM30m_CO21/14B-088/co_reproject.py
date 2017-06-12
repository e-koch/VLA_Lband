
from astropy.io import fits
from reproject import reproject_interp
from spectral_cube.lower_dimensional_structures import Projection
from astropy.wcs import WCS

from analysis.paths import fourteenB_HI_data_path, iram_co21_data_path
from analysis.constants import moment0_name

'''
Reproject the CO to match the HI.

Save one version with the CO smoothed to match the HI

'''

hi_header = fits.getheader(fourteenB_HI_data_path(moment0_name))

co_ico = fits.open(iram_co21_data_path("m33.ico.fits"))[0]

proj = Projection(co_ico.data.squeeze(),
                  wcs=WCS(co_ico.header).dropaxis(3).dropaxis(2))

rep_array, footprint = reproject_interp(proj.hdu, hi_header)

new_header = hi_header.copy()

# Need to change some of the properties before saving
keys = ["BMAJ", "BMIN", "BPA", "BUNIT"]
for key in keys:
    new_header[key] = co_ico.header[key]

new_hdu = fits.PrimaryHDU(rep_array, header=new_header)
new_hdu.writeto(iram_co21_data_path("m33.ico.hireprojection.fits",
                                    no_check=True), clobber=True)
