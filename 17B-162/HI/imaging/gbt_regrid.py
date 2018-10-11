
'''
Regrid the GBT data to match the VLA HI data.
'''

from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar
import numpy as np
import os

from cube_analysis.io_utils import create_huge_fits

from paths import (seventeenB_HI_data_02kms_path,
                   seventeenB_HI_data_1kms_path, data_path)

run_02kms = True
run_1kms = False

gbt_path = os.path.join(data_path, "GBT")

if not os.path.exists(os.path.join(gbt_path, '17B-162_items')):
    os.mkdir(os.path.join(gbt_path, '17B-162_items'))

# Ta* to T_mb as per @low-sky
Tmb_conv = 1.052

cube = SpectralCube.read(os.path.join(gbt_path, "m33_gbt_vlsr_highres.fits"))

if run_02kms:
    # Load the non-pb masked cube
    vla_cube = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.fits"))

    save_name = os.path.join(gbt_path, "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_02kms_spectralregrid.fits")

    # Spectral interpolation, followed by reprojection.
    if not os.path.exists(save_name):

        cube = cube.spectral_interpolate(vla_cube.spectral_axis)

        if cube._is_huge:
            output_fits = create_huge_fits(save_name, cube.header,
                                           return_hdu=True)

            for chan in ProgressBar(cube.shape[0]):
                output_fits[0].data[chan] = cube[chan].value * Tmb_conv
            output_fits.flush()
            output_fits.close()
        else:
            (cube * Tmb_conv).write(save_name, overwrite=True)
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
    out_path = os.path.join(gbt_path, '17B-162_items')

    save_name = os.path.join(gbt_path,
                             "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_02kms.fits")
    output_fits = create_huge_fits(save_name, new_header, return_hdu=True)

    targ_header = vla_cube[0].header
    for chan in ProgressBar(cube.shape[0]):
        reproj_chan = \
            cube[chan].reproject(targ_header).value.astype(np.float32)
        output_fits[0].data[chan] = reproj_chan
        if chan % 100 == 0:
            output_fits.flush()
    output_fits.close()

if run_1kms:
    # Load the non-pb masked cube
    vla_cube = SpectralCube.read(seventeenB_HI_data_1kms_path("M33_14_17B_HI_contsub_width_1kms.pbcor.image.fits"))

    save_name = os.path.join(gbt_path,
                             "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms_spectralregrid.fits")

    # Spectral interpolation, followed by reprojection.
    if not os.path.exists(save_name):

        cube = cube.spectral_interpolate(vla_cube.spectral_axis)

        if cube._is_huge:
            output_fits = create_huge_fits(save_name, cube.header,
                                           return_hdu=True)

            for chan in ProgressBar(cube.shape[0]):
                output_fits[0].data[chan] = cube[chan].value * Tmb_conv
            output_fits.flush()
            output_fits.close()
        else:
            (cube * Tmb_conv).write(save_name, overwrite=True)
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
    out_path = os.path.join(gbt_path, '17B-162_items')

    save_name = os.path.join(gbt_path,
                             "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms.fits")
    output_fits = create_huge_fits(save_name, new_header, return_hdu=True)

    targ_header = vla_cube[0].header
    for chan in ProgressBar(cube.shape[0]):
        reproj_chan = \
            cube[chan].reproject(targ_header).value.astype(np.float32)
        output_fits[0].data[chan] = reproj_chan
        if chan % 100 == 0:
            output_fits.flush()
    output_fits.close()
