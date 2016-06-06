
from astropy.io import fits
from reproject import reproject_interp
from spectral_cube.lower_dimensional_structures import Projection
import os
from astropy.wcs import WCS

'''
Reproject the CO to match the HI.

Not taking the beam difference into account here since they are quite close.
'''

hi_header = fits.getheader("/home/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits")

co_path = "/media/eric/Data_3/M33/co21"

co_ico = fits.open(os.path.join(co_path, "m33.ico.fits"))[0]

proj = Projection(co_ico.data.squeeze(), wcs=WCS(co_ico.header).dropaxis(3).dropaxis(2))

rep_array, footprint = reproject_interp(proj.hdu, hi_header)

new_header = hi_header.copy()

# Need to change some of the properties before saving
keys = ["BMAJ", "BMIN", "BPA", "BUNIT"]
for key in keys:
    new_header[key] = co_ico.header[key]

new_hdu = fits.PrimaryHDU(rep_array, header=new_header)
new_hdu.writeto(os.path.join(co_path, "m33.ico.hireprojection.fits"))
