
'''
Regrid the GBT data to match the VLA HI data.
'''

from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar
import numpy as np
import astropy.io.fits as fits
import os

from paths import fourteenB_HI_data_path, data_path


# Load the non-pb masked cube
vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))

gbt_path = os.path.join(data_path, "GBT")

cube = SpectralCube.read(os.path.join(gbt_path, "m33_gbt_vlsr_highres.fits"))

# Spectral interpolation, followed by reprojection.
cube = cube.spectral_interpolate(vla_cube.spectral_axis)
save_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid.fits")
# Cut to the VLA map extents. Technically need to reproject, but this by-eye
# pixel selection is just slightly beyond the VLA data
# cube[:, 20:80, 40:80].write(save_name)
cube.write(save_name)

# Make the reprojected header
new_header = cube.header.copy()
new_header["NAXIS"] = 3
new_header["NAXIS1"] = vla_cube.shape[2]
new_header["NAXIS2"] = vla_cube.shape[1]
new_header["NAXIS3"] = vla_cube.shape[0]
kwarg_skip = ['TELESCOP', 'BUNIT', 'INSTRUME']
for key in cube.header:
    if key == 'HISTORY':
        continue
    if key in vla_cube.header:
        if "NAXIS" in key:
            continue
        if key in kwarg_skip:
            continue
        new_header[key] = vla_cube.header[key]
new_header.update(cube.beam.to_header_keywords())
new_header["BITPIX"] = -32
# We're going to convert to Tmb below
new_header.comments['BUNIT'] = 'Tmb'

# Build up the reprojected cube per channel
save_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088.fits")
output_fits = fits.StreamingHDU(save_name, new_header)

fill_plane = np.zeros_like(vla_cube[0].value, dtype=np.float32) * np.NaN
for chan in ProgressBar(cube.shape[0]):
    output_fits.write(fill_plane)
output_fits.close()

output_fits = fits.open(save_name, mode='update')

targ_header = vla_cube[0].header
for chan in ProgressBar(cube.shape[0]):
    reproj_chan = \
        cube[chan].reproject(targ_header).value.astype(np.float32)
    # Ta* to T_mb as per @low-sky
    output_fits[0].data[chan] = reproj_chan * 1.052
    if chan % 20 == 0:
        output_fits.flush()
output_fits.close()
