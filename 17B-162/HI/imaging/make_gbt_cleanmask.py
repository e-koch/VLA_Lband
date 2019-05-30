
'''
Using really large scales in the multiscale clean makes is susceptible to
diverging. Set the large-scale mask based on the GBT data multiplied by
the pb mask.

Use threshold mask produced by gbt_regrid.py
'''

import os
import numpy as np
# from spectral_cube import SpectralCube
# from spectral_cube.io.casa_masks import make_casa_mask

from tasks import importfits

from paths import data_path

run_02kms = False
run_1kms = True

gbt_path = os.path.join(data_path, "GBT")

mask_save_name = os.path.join(gbt_path,
                              "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms_mask.fits")
mask_out_name = os.path.join(gbt_path,
                             "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms_mask.image")

importfits(fitsimage=mask_save_name, imagename=mask_out_name, defaultaxes=True,
           defaultaxesvalues=["", "", "", "I"])

# cube = SpectralCube.read(mask_save_name)
# make_casa_mask(cube, mask_out_name,
#                append_to_image=False)
# del cube

# Split off each channel
out_path = os.path.join(gbt_path, "17B-162_items/mask_1kms_chans/")

if not os.path.exists(out_path):
    os.mkdir(out_path)

# Open image.
ia.open(mask_out_name)

# Find the spectral axis
csys = ia.coordsys()

# Check this name??
specaxis_name = "Frequency"
spec_axis = np.where(np.asarray(csys.names()) == specaxis_name)[0][0]

cube_shape = list(ia.shape())
ndims = len(cube_shape)
nchan = cube_shape[spec_axis]

lower_corner = [0] * ndims
upper_corner = cube_shape

for chan in range(nchan):

    casalog.post("Channel {}".format(chan))

    # Set the channel
    lower_corner[spec_axis] = chan
    upper_corner[spec_axis] = chan

    box = rg.box(lower_corner, upper_corner)

    # Now make sliced image
    chan_name = "{0}_channel_{1}".format(mask_out_name.split("/")[-1], chan)
    im_slice = ia.subimage(os.path.join(out_path, chan_name),
                           box)
    im_slice.done()

ia.close()
