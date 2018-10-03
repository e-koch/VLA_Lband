
'''
Combine individually-imaged channels into a cube.

Run in CASA.
'''

import sys
from glob import glob

from taskinit import iatool

ia = iatool()

path_to_data = sys.argv[-4]

filename = sys.argv[-3]

# Check for the expected number of images
num_imgs = int(sys.argv[-2])

suffix = sys.argv[-1]

suffixes = ['mask', 'model', 'pb', 'psf', 'residual', 'residual_init', 'image',
            'image.pbcor', 'sumwt', 'weight']

if suffix not in suffixes:
    raise NameError("suffix {0} is not a valid output file type from tclean.")

casalog.post("Assembling {} cube".format(suffix))

images = []

for num in range(num_imgs):

    chan_img = glob("{0}/channel_{1}/{2}_channel*.{3}"
                    .format(path_to_data, num, filename, suffix))
    if len(chan_img) == 1:
        chan_img.append(chan_img[0])

if len(images) == 0:
    casalog.post("No images found for {}".format(suffix))
    continue

if len(images) != num_imgs:
    casalog.post("Number of images found ({0}) does not match"
                 " expected number ({1}) for {2}. Skipping cube creation."
                 .format(len(images), num_imgs, suffix))
    continue

cubename = "{0}/{1}.{2}".format(path_to_data, filename, suffix)

ia.imageconcat(outfile=cubename, infiles=images, reorder=False,
               overwrite=True)
ia.close()
