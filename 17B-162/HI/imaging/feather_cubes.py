
'''
Create feathered VLA cubes with the 3 different SD datasets.
'''

from spectral_cube import SpectralCube
import os
from os.path import join as osjoin
from astropy import log
import scipy.ndimage as nd
import numpy as np

from cube_analysis.feather_cubes import feather_cube

from paths import (seventeenB_HI_data_02kms_path,
                   seventeenB_HI_data_1kms_path,
                   data_path)
from constants import hi_freq

# Set which of the cubes to feather
run_gbt_02kms = True
run_gbt_1kms = False

num_cores = 4
chunk = 8


def taper_weights(mask, sigma, nsig_cut=3):
    '''
    This needs to be moved to uvcombine.
    '''

    dist = nd.distance_transform_edt(mask)

    gauss_dists = np.where(np.logical_and(dist < nsig_cut * sigma, dist > 0.))
    flat_dists = np.where(dist >= nsig_cut * sigma)

    weight_arr = np.zeros_like(mask, dtype=float)

    weight_arr[gauss_dists] = \
        np.exp(- (dist[gauss_dists] - nsig_cut * sigma)**2 / (2 * sigma**2))
    weight_arr[flat_dists] = 1.

    return weight_arr


if run_gbt_02kms:
    log.info("Feathering with 0.2 km/s GBT")

    # Load the non-pb masked cube
    vla_cube = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.image.pbcor.fits"))

    pb_cube = SpectralCube.read(seventeenB_HI_data_02kms_path("M33_14B_17B_HI_contsub_width_02kms.pb.fits"))
    # PB minimally changes over the frequency range. So just grab one plane
    pb_plane = pb_cube[0]

    # Smoothly taper data at the mosaic edge. This weight array tapers to
    # exp(-5^2/2)~4e-6 at the pb cut-off of 0.2.
    weight = taper_weights(np.isfinite(pb_plane), 30, nsig_cut=5)

    gbt_path = osjoin(data_path, "GBT")
    gbt_name = osjoin(gbt_path, "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_02kms_spectralregrid.fits")
    gbt_cube = SpectralCube.read(gbt_name)

    output_path = osjoin(data_path,
                         "VLA/17B-162/HI/full_imaging_02kms_wGBT/")
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    save_name = osjoin(output_path,
                       "M33_14B_17B_HI_contsub_width_02kms.image.pbcor.GBT_feathered.fits")

    feather_cube(vla_cube, gbt_cube, restfreq=hi_freq, save_feather=True,
                 save_name=save_name, num_cores=num_cores,
                 weights=weight, chunk=chunk, verbose=False)

if run_gbt_1kms:
    log.info("Feathering with 1 km/s GBT")

    # Load the non-pb masked cube
    vla_cube = SpectralCube.read(seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.image.pbcor.fits"))

    pb_cube = SpectralCube.read(seventeenB_HI_data_1kms_path("M33_14B_17B_HI_contsub_width_1kms.pb.fits"))
    # PB minimally changes over the frequency range. So just grab one plane
    pb_plane = pb_cube[0]

    # Smoothly taper data at the mosaic edge. This weight array tapers to
    # exp(-5^2/2)~4e-6 at the pb cut-off of 0.2.
    weight = taper_weights(np.isfinite(pb_plane), 30, nsig_cut=5)

    gbt_path = osjoin(data_path, "GBT")
    gbt_name = osjoin(gbt_path,
                      "17B-162_items/m33_gbt_vlsr_highres_Tmb_17B162_1kms.fits")
    gbt_cube = SpectralCube.read(gbt_name)

    output_path = osjoin(data_path, "VLA/17B-162/HI/full_imaging_1kms_wGBT/")
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    save_name = osjoin(output_path,
                       "M33_14B_17B_HI_contsub_width_1kms.image.pbcor.GBT_feathered.fits")

    feather_cube(vla_cube, gbt_cube, restfreq=hi_freq, save_feather=True,
                 save_name=save_name, num_cores=num_cores,
                 weights=weight, chunk=chunk, verbose=False)
