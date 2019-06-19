
'''
Using really large scales in the multiscale clean makes is susceptible to
diverging. Set the large-scale mask based on the GBT data multiplied by
the pb mask.

Use threshold mask produced by gbt_regrid.py
'''

import os
import numpy as np
from spectral_cube import SpectralCube
from reproject import reproject_interp
from radio_beam import Beam
# from spectral_cube.io.casa_masks import make_casa_mask
import astropy.units as u
from astropy.io import fits
from scipy import ndimage as nd
from astropy.wcs import WCS
from astropy.utils.console import ProgressBar

from cube_analysis.io_utils import create_huge_fits

from paths import (data_path, fourteenB_wGBT_HI_file_dict,
                   fourteenB_HI_data_wGBT_path,
                   seventeenB_HI_data_1kms_path,
                   )

run_02kms = False
run_1kms = True


def vel_to_freq(vel_or_freq, rest_freq=1.42040575177 * u.GHz,
                unit=u.Hz):
    '''
    Using radio velocity here.
    '''
    equiv = u.doppler_radio(rest_freq)

    return vel_or_freq.to(unit, equiv)


if run_02kms:
    raise NotImplementedError("")

if run_1kms:

    # Try using the 14B source mask as a guide for the clean mask.

    mask_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["Source_Mask"])

    # Good enough for the mask
    beam = Beam(19 * u.arcsec)
    mask_cube = mask_cube.with_beam(beam)

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

    save_name = fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.pbcov_gt_0.5_masked_source_mask_17B_1kms_spectral_match.fits",
                                            no_check=True)

    if not os.path.exists(save_name):

        mask_cube = mask_cube.spectral_interpolate(vel_axis)

        # if cube._is_huge:
        #     output_fits = create_huge_fits(save_name, mask_cube.header,
        #                                    return_hdu=True)

        #     for chan in ProgressBar(cube.shape[0]):
        #         output_fits[0].data[chan] = cube[chan].value * Tmb_conv
        #     output_fits.flush()
        #     output_fits.close()
        # else:
        mask_cube.write(save_name, overwrite=True)
    else:
        mask_cube = SpectralCube.read(save_name)
        mask_cube = mask_cube.with_beam(beam)

    # Smooth out the mask and fill in small holes
    beam_element = Beam(19 * 10 * u.arcsec).as_tophat_kernel((3 * u.arcsec).to(u.deg)).array > 0

    new_header = mask_cube.header.copy()
    new_header["NAXIS"] = 3
    new_header["NAXIS1"] = vla_pbmask.shape[1]
    new_header["NAXIS2"] = vla_pbmask.shape[0]
    new_header["NAXIS3"] = nchan
    new_header['CRVAL3'] = vel_axis[0].value
    kwarg_skip = ['TELESCOP', 'BUNIT', 'INSTRUME']
    for key in mask_cube.header:
        if key == 'HISTORY' or key == "COMMENT":
            continue
        if key in vla_spat_hdr:
            if "NAXIS" in key:
                continue
            if key in kwarg_skip:
                continue
            new_header[key] = vla_spat_hdr[key]
    new_header.update(mask_cube.beam.to_header_keywords())

    new_header['BITPIX'] = 8
    new_header['BUNIT'] = 'bool'

    save_name = "{}_17B_reproj.fits".format(fourteenB_wGBT_HI_file_dict["Source_Mask"].rstrip(".fits"))
    output_fits = create_huge_fits(save_name, new_header, return_hdu=True)

    targ_header = WCS(vla_spat_hdr).celestial.to_header()
    targ_header["NAXIS"] = 2
    targ_header["NAXIS1"] = vla_pbmask.shape[1]
    targ_header["NAXIS2"] = vla_pbmask.shape[0]

    for chan in ProgressBar(nchan):
        mask_plane = mask_cube[chan].value.astype(bool)

        mask_plane = nd.binary_dilation(mask_plane, beam_element)

        # Now reproject
        mask_plane = reproject_interp((mask_plane.astype(float), mask_cube[chan].header), targ_header)[0] > 1e-3

        # Apply the pbmask
        mask_plane = np.logical_and(mask_plane, vla_pbmask)

        output_fits[0].data[chan] = mask_plane.astype(">i2")
        if chan % 50 == 0:
            output_fits.flush()

    output_fits.close()

    # gbt_path = os.path.join(data_path, "GBT")

    # mask_save_name = os.path.join(gbt_path,
    #                               "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms_mask.fits")
    # mask_out_name = os.path.join(gbt_path,
    #                              "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms_mask.image")
