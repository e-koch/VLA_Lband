
'''
Combine individually-imaged channels into a cube.

Run in CASA.
'''

import sys
from glob import glob
import os

from taskinit import iatool

ia = iatool()

path_to_data = sys.argv[-3]

filename = sys.argv[-2]

# Check for the expected number of images
num_imgs = int(sys.argv[-1])

suffixes = ['mask', 'model', 'pb', 'psf', 'residual', 'residual_init', 'image',
            'image.pbcor', 'sumwt', 'weight']

for suff in suffixes:

    casalog.post("Assembling {} cube".format(suff))

    images = []

    for num in range(num_imgs):

        chan_img = glob("{0}/{1}_channel*.{2}".format(path_to_data, filename,
                                                      suff))
        if chan_img == 1:
            chan_img.append(chan_img[0])

    if len(images) == 0:
        casalog.post("No images found for {}".format(suff))
        continue

    if len(images) != num_imgs:
        casalog.post("Number of images found ({0}) does not match"
                     " expected number ({1}) for {2}. Skipping cube creation."
                     .format(len(images), num_imgs, suff))
        continue

    cubename = "{0}/{1}.{2}".format(path_to_data, filename, suff)

    ia.imageconcat(outfile=cubename, infiles=images, reorder=True,
                   overwrite=True)
    ia.close()
