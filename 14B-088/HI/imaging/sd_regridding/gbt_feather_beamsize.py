
'''
How much is the feathered image affected by assuming a wrong beam size?
'''

from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import os
from corner import hist2d
from radio_beam import Beam
import astropy.units as u
from astropy.visualization import hist
from uvcombine.uvcombine import feather_simple, feather_compare
import numpy as np
import scipy.ndimage as nd

from paths import fourteenB_HI_data_path, data_path
from constants import hi_freq

vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))

gbt_path = os.path.join(data_path, "GBT")
gbt_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid.fits")
gbt_cube = SpectralCube.read(gbt_name)

gbt_fullregrid_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088.fits")
gbt_fullregrid_cube = SpectralCube.read(gbt_fullregrid_name)

beam_fwhm = lambda diam: ((1.2 * 21 * u.cm) / diam.to(u.cm)) * u.rad

chan = 500

gbt_plane = gbt_cube[chan]
gbt_fr_plane = gbt_fullregrid_cube[chan]
vla_plane = vla_cube[chan].to(u.K, vla_cube.beams[chan].jtok_equiv(hi_freq))

feather_80 = feather_simple(vla_plane.hdu, gbt_plane.hdu,
                            lowresfwhm=beam_fwhm(80 * u.m).to(u.arcsec))

feather_90 = feather_simple(vla_plane.hdu, gbt_plane.hdu,
                            lowresfwhm=beam_fwhm(90 * u.m).to(u.arcsec))

feather_90_fr = feather_simple(vla_plane.hdu, gbt_fr_plane.hdu,
                               lowresfwhm=beam_fwhm(90 * u.m).to(u.arcsec))

feather_100 = feather_simple(vla_plane.hdu, gbt_plane.hdu,
                             lowresfwhm=beam_fwhm(100 * u.m).to(u.arcsec))

mask = gbt_fr_plane.value > 2

vla_beam_kernel = vla_plane.beam.as_tophat_kernel(vla_plane.header["CDELT2"]).array > 0
vla_mask = np.isfinite(vla_plane)
vla_mask = nd.binary_erosion(vla_mask, vla_beam_kernel, iterations=10)

plt.plot([feather_80.real[vla_mask].sum() / gbt_fr_plane[vla_mask].sum().value,
          feather_90.real[vla_mask].sum() / gbt_fr_plane[vla_mask].sum().value,
          feather_100.real[vla_mask].sum() / gbt_fr_plane[vla_mask].sum().value])
