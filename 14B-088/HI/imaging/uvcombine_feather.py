
'''
Use uvcombine to feather with SD data
'''

from uvcombine import feather_simple
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import beams_to_bintable
from radio_beam import Beam
from astropy.utils.console import ProgressBar
import numpy as np
import astropy.io.fits as fits
import astropy.units as u
import os

from paths import fourteenB_HI_data_path, data_path
from constants import hi_freq


# Load the non-pb masked cube
vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))

gbt_path = os.path.join(data_path, "GBT")
gbt_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid.fits")
gbt_cube = SpectralCube.read(gbt_name)

output_path = os.path.join(data_path, "VLA/14B-088/HI/full_imaging_wGBT/")

save_name = os.path.join(output_path, "M33_14B-088_HI.clean.image.GBT_highres_feathered.fits")
output_fits = fits.StreamingHDU(save_name, vla_cube.header)

fill_plane = np.zeros_like(vla_cube[0].value, dtype=np.float32) * np.NaN
for chan in ProgressBar(vla_cube.shape[0]):
    output_fits.write(fill_plane)
output_fits.close()

output_fits = fits.open(save_name, mode='update')

wave = hi_freq.to(u.m, equivalencies=u.spectral())
effdish_diam = 80 * u.m
targ_beam = Beam(1.22 * (wave / effdish_diam) * u.rad)

vla_mask = np.isfinite(vla_cube[0])

for chan in ProgressBar(vla_cube.shape[0]):
    vla_plane = vla_cube[chan].hdu
    vla_plane.header['BUNIT'] = "Jy/beam"
    gbt_plane = gbt_cube[-chan].to(u.Jy, equivalencies=targ_beam.jtok_equiv(hi_freq)).hdu
    gbt_plane.header['BUNIT'] = "Jy/beam"
    gbt_plane.header.update(targ_beam.to_header_keywords())
    feathered = feather_simple(vla_plane, gbt_plane).real.astype(np.float32)
    feathered[~vla_mask] = np.NaN
    output_fits[0].data[chan] = feathered
    if chan % 20 == 0:
        output_fits.flush()

# Finally, append the beam table
output_fits.append(beams_to_bintable(vla_cube.beams))

output_fits.close()
