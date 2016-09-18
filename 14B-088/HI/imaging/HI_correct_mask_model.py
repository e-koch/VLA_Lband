
'''
\Swap the spectral and stokes axes. Needed due to issue in regridding function
'''

import sys
from astropy.io import fits


hdu = fits.open(sys.argv[1], mode='update')

hdu[0].data = hdu[0].data.swapaxes(0, 1)

execfile("/home/eric/Dropbox/code_development/ewky_scripts/header_swap_axis.py")

hdu[0].header = header_swapaxes(hdu[0].header, 2, 3)

hdu.flush()

hdu.close()
