
'''
Dilate the existing clean mask.
'''

from astropy.io import fits
from skimage.morphology import disk
from scipy import ndimage as nd


mask = fits.open("M33_14B-088_HI_mask_modified.fits", mode='update')
pbcov = fits.getdata("M33_14B-088_pbcor.fits")

pb_thresh = 0.2
pb_mask = pbcov > pb_thresh

struct = disk(100)

for i in xrange(1231):
    print(i)
    mask[0].data[i, 0, :, :] = \
        nd.binary_dilation(mask[0].data[i, 0, :, :], structure=struct)
    mask[0].data[i, 0, :, :] *= pb_mask.squeeze()

    mask.flush()

mask.close()
