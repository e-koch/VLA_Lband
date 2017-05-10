
'''
Regrid the EBHIS data to match the VLA HI data.
'''

from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar
import numpy as np
import os

from cube_analysis.io_utils import create_huge_fits

from paths import fourteenB_HI_data_path, data_path


# Load the non-pb masked cube
vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))

ebhis_path = os.path.join(data_path, "EBHIS")

# Spectral interpolation, followed by reprojection.
save_name = os.path.join(ebhis_path, "14B-088_items/m33_ebhis_14B088_spectralregrid.fits")

# Only load in the original data if the spectrally-interpolated cube hasn't
# been saved
if not os.path.exists(save_name):

    cube = SpectralCube.read(os.path.join(ebhis_path, "CAR_C02.fit"))

    # Cut to the spatial extents, since this chunk of EBHIS covers so much more
    cube = cube.subcube(xlo=vla_cube.longitude_extrema[1],
                        xhi=vla_cube.longitude_extrema[0],
                        ylo=vla_cube.latitude_extrema[0],
                        yhi=vla_cube.latitude_extrema[1])

    cube = cube.spectral_interpolate(vla_cube.spectral_axis)
    # Write out the spectrally interpolated cube
    cube.write(save_name)
else:
    cube = SpectralCube.read(save_name)

# Make the reprojected header
new_header = cube.header.copy()
new_header["NAXIS"] = 3
new_header["NAXIS1"] = vla_cube.shape[2]
new_header["NAXIS2"] = vla_cube.shape[1]
new_header["NAXIS3"] = vla_cube.shape[0]
kwarg_skip = ['TELESCOP', 'BUNIT', 'INSTRUME']
for key in cube.header:
    if key == 'HISTORY' or key == 'COMMENT':
        continue
    if key in vla_cube.header:
        if "NAXIS" in key:
            continue
        if key in kwarg_skip:
            continue
        new_header[key] = vla_cube.header[key]
new_header.update(cube.beam.to_header_keywords())
new_header["BITPIX"] = -32

# Build up the reprojected cube per channel
save_name = os.path.join(ebhis_path, "14B-088_items/m33_ebhis_14B088.fits")
output_fits = create_huge_fits(save_name, new_header, return_hdu=True)

targ_header = vla_cube[0].header
for chan in ProgressBar(cube.shape[0]):
    proj = cube[chan].reproject(targ_header).value.astype(np.float32)
    output_fits[0].data[chan] = proj
    if chan % 100 == 0:
        output_fits.flush()
output_fits.close()
