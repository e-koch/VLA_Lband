
'''
Combine individually-imaged channels into a cube.

Run in CASA.
'''

import sys
from glob import glob
import socket

from taskinit import iatool

ia = iatool()

# https://github.com/e-koch/VLA_Lband/blob/master/CASA_functions/imaging_utils.py
if socket.gethostname().lower() == 'segfault':
    execfile("/home/ekoch/ownCloud/code_development/VLA_Lband/CASA_functions/imaging_utils.py")
else:
    execfile("/home/ekoch/code/VLA_Lband/CASA_functions/imaging_utils.py")

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

    chan_img = glob("{0}/{1}_channel_{2}.{3}"
                    .format(path_to_data, filename, num, suffix))
    if len(chan_img) == 1:
        images.append(chan_img[0])

if len(images) == 0:
    casalog.post("No images found for {}".format(suffix))
    sys.exit(1)

if len(images) != num_imgs:
    casalog.post("Number of images found ({0}) does not match"
                 " expected number ({1}) for {2}. Skipping cube creation."
                 .format(len(images), num_imgs, suffix))
    sys.exit(1)

cubename = "{0}/{1}.{2}".format(path_to_data, filename, suffix)

append_to_cube(path_to_data, "{}_channel".format(filename),
               suffix, num_imgs,
               cubename, chunk_size=250,
               delete_chunk_cubes=True,
               concat_kwargs={'relax': True, 'reorder': False,
                              'overwrite': True})

casalog.post("Look! I made a {} cube!".format(suffix))
