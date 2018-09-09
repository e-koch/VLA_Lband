
'''
What are the differences between clean and tclean?

Run tclean_test.py first.
'''

from spectral_cube import Projection
from astropy.io import fits


clean_dirty = Projection.from_hdu(fits.open("clean_test_dirty.image.fits")[0])
clean_img = Projection.from_hdu(fits.open("clean_test.image.fits")[0])
clean_model = Projection.from_hdu(fits.open("clean_test.model.fits")[0])
clean_residual = Projection.from_hdu(fits.open("clean_test.residual.fits")[0])

tclean_dirty = Projection.from_hdu(fits.open("tclean_test_dirty.image.fits")[0])
tclean_img = Projection.from_hdu(fits.open("tclean_test.image.fits")[0])
tclean_model = Projection.from_hdu(fits.open("tclean_test.model.fits")[0])
tclean_residual = Projection.from_hdu(fits.open("tclean_test.residual.fits")[0])

tclean_casa5_img = Projection.from_hdu(fits.open("tclean_test_casa5.0.image.fits")[0])
tclean_casa5_model = Projection.from_hdu(fits.open("tclean_test_casa5.0.model.fits")[0])
tclean_casa5_residual = Projection.from_hdu(fits.open("tclean_test_casa5.0.residual.fits")[0])


# So, what's different or the same?
