
'''
Regrid the GBT data to match the VLA HI data.

Generated from spatial info + a constructed spectral axis.
Also creates a clean mask based on a visually-determined
threshold.
'''

from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import SpectralCube
from radio_beam import Beam
from astropy.utils.console import ProgressBar
import numpy as np
import os
import astropy.units as u

from cube_analysis.io_utils import create_huge_fits

from paths import (seventeenB_HI_data_02kms_path,
                   seventeenB_HI_data_1kms_path, data_path)
from constants import hi_freq

run_02kms = False
run_1kms = True

gbt_path = os.path.join(data_path, "GBT")

if not os.path.exists(os.path.join(gbt_path, '17B-162_items')):
    os.mkdir(os.path.join(gbt_path, '17B-162_items'))

# Ta* to T_mb as per @low-sky
Tmb_conv = 1.052

cube = SpectralCube.read(os.path.join(gbt_path, "m33_gbt_vlsr_highres.fits"))

# Update the beams to be the optimal found for the 14B data
beam_fwhm = lambda diam: ((1.18 * hi_freq.to(u.cm, u.spectral())) /
                           diam.to(u.cm)) * u.rad
gbt_eff_beam = beam_fwhm(87.5 * u.m).to(u.deg)


# Start and stop velocities (as used in match_and_split.py)

# in km/s
start_vel = -330
end_vel = -50


def vel_to_freq(vel_or_freq, rest_freq=1.42040575177 * u.GHz,
                unit=u.Hz):
    '''
    Using radio velocity here.
    '''
    equiv = u.doppler_radio(rest_freq)

    return vel_or_freq.to(unit, equiv)


if run_02kms:
    # Load the non-pb masked cube
    # vla_cube = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.fits"))

    # Load the non-pb masked cube
    # vla_cube = SpectralCube.read(seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.image.pbcor.fits"))
    vla_pbmask = fits.getdata(seventeenB_HI_data_1kms_path("17B_14B_spatial_pbmask.fits")) > 0

    vla_spat_hdr = fits.Header.fromtextfile(seventeenB_HI_data_1kms_path("17B_14B_spatial_header.txt"))

    print("YOU NEED TO UPDATE THIS PART")
    print(argh)
    # Hard code in properties to make the spectral axis
    chan_width = 0.21

    # vel_0 = XXX
    # del_vel = np.abs(vla_spat_hdr['CDELT3'])
    # nchan = 1359

    save_name = os.path.join(gbt_path, "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_02kms_spectralregrid.fits")

    # Spectral interpolation, followed by reprojection.
    if not os.path.exists(save_name):

        cube = cube.spectral_interpolate(vla_cube.spectral_axis)
        cube = cube.with_beam(Beam(gbt_eff_beam))

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
    # vla_cube = SpectralCube.read(seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.image.pbcor.fits"))
    vla_pbmask = fits.getdata(seventeenB_HI_data_1kms_path("17B_14B_spatial_pbmask.fits")) > 0

    vla_spat_hdr = fits.Header.fromtextfile(seventeenB_HI_data_1kms_path("17B_14B_spatial_header.txt"))

    # Hard code in properties to make the spectral axis

    # Define axis in frequency. Then convert to V_rad
    freq_0 = 1420642320.9823689 * u.Hz
    # This cube averages over 5 channels
    del_freq = 976.5571551322937 * 5 * u.Hz
    nchan = 272

    freq_axis = np.arange(nchan) * del_freq + freq_0

    vel_axis = vel_to_freq(freq_axis, unit=u.m / u.s)

    save_name = os.path.join(gbt_path,
                             "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms_spectralregrid.fits")

    # Spectral interpolation, followed by reprojection.
    if not os.path.exists(save_name):

        cube = cube.spectral_interpolate(vel_axis)
        cube = cube.with_beam(Beam(gbt_eff_beam))

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
    new_header["NAXIS1"] = vla_pbmask.shape[1]
    new_header["NAXIS2"] = vla_pbmask.shape[0]
    new_header["NAXIS3"] = nchan
    new_header['CRVAL3'] = vel_axis[0].value
    kwarg_skip = ['TELESCOP', 'BUNIT', 'INSTRUME']
    for key in cube.header:
        if key == 'HISTORY':
            continue
        if key in vla_spat_hdr:
            if "NAXIS" in key:
                continue
            if key in kwarg_skip:
                continue
            new_header[key] = vla_spat_hdr[key]
    new_header.update(cube.beam.to_header_keywords())
    new_header["BITPIX"] = -32
    # We're going to convert to Tmb below
    new_header.comments['BUNIT'] = 'Tmb'

    # Build up the reprojected cube per channel
    out_path = os.path.join(gbt_path, '17B-162_items')

    save_name = os.path.join(gbt_path,
                             "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms.fits")
    output_fits = create_huge_fits(save_name, new_header, return_hdu=True)

    targ_header = WCS(vla_spat_hdr).celestial.to_header()
    targ_header["NAXIS"] = 2
    targ_header["NAXIS1"] = vla_pbmask.shape[1]
    targ_header["NAXIS2"] = vla_pbmask.shape[0]

    for chan in ProgressBar(cube.shape[0]):
        reproj_chan = \
            cube[chan].reproject(targ_header).value.astype(np.float32)
        output_fits[0].data[chan] = reproj_chan
        if chan % 100 == 0:
            output_fits.flush()
    output_fits.close()

    del output_fits
    del cube

    # Also create a clean mask version from the regridded data.

    thresh = 1.0 * u.K

    # Re-open the output an apply the brightness threshold to create a mask
    new_header['BITPIX'] = 8
    new_header['BUNIT'] = 'bool'
    mask_save_name = os.path.join(gbt_path,
                             "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms_mask.fits")
    output_fits = create_huge_fits(mask_save_name, new_header, return_hdu=True)

    hdu = fits.open(save_name, mode='denywrite')

    for chan in ProgressBar(nchan):
        mask_chan = hdu[0].data[chan] > thresh.value
        # Apply the pbmask
        mask_chan = np.logical_and(mask_chan, vla_pbmask)

        output_fits[0].data[chan] = mask_chan.astype(">i2")
        if chan % 50 == 0:
            output_fits.flush()

    output_fits.close()
