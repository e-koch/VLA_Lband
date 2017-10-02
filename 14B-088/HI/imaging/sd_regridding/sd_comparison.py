
'''
Compare the regridded versions of the SD datasets.
'''

from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import os
from corner import hist2d
from radio_beam import Beam
import astropy.units as u
import numpy as np

from paths import fourteenB_HI_data_path, data_path
from galaxy_params import gal

# Load in the 4 cubes and run.

vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))


arecibo_path = os.path.join(data_path, "Arecibo")

# Spectral interpolation, followed by reprojection.
arecibo_name = \
    os.path.join(arecibo_path,
                 "14B-088_items_new/m33_arecibo_14B088.fits")
arecibo_cube = SpectralCube.read(arecibo_name)

ebhis_path = os.path.join(data_path, "EBHIS")

# Spectral interpolation, followed by reprojection.
ebhis_name = os.path.join(ebhis_path, "14B-088_items/m33_ebhis_14B088.fits")
ebhis_cube = SpectralCube.read(ebhis_name)


gbt_path = os.path.join(data_path, "GBT")
gbt_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088.fits")
gbt_cube = SpectralCube.read(gbt_name)

gbt_lowres_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_Tmb_14B088.fits")
gbt_lowres_cube = SpectralCube.read(gbt_lowres_name)

# Compare total emission in the cubes.
vla_mask = np.isfinite(vla_cube[0])

arecibo_sum = arecibo_cube.with_mask(vla_mask).sum()
ebhis_sum = ebhis_cube.with_mask(vla_mask).sum()
gbt_sum = gbt_cube.with_mask(vla_mask).sum()
gbt_lowres_sum = gbt_lowres_cube.with_mask(vla_mask).sum()

plt.plot(arecibo_sum, ebhis_sum, gbt_sum, gbt_lowres_sum)

# Compare intensities in one plane

# arecibo_plane = arecibo_cube[500]
# ebhis_plane = ebhis_cube[500]
# gbt_plane = gbt_cube[500]
# gbt_plane[np.isnan(gbt_plane)] = 0.0 * u.K
# gbt_lowres_plane = gbt_lowres_cube[500]

# # Convolve GBT to match EBHIS
# beam_fwhm = lambda diam: ((1.2 * 21 * u.cm) / diam.to(u.cm)) * u.rad
# gbt_90m_beam = Beam(beam_fwhm(90 * u.m))

# gbt_plane._beam = gbt_90m_beam

# gbt_plane_convolved = gbt_plane.convolve_to(ebhis_plane.beam)


# gbt_100m_beam = Beam(beam_fwhm(100 * u.m))
# gbt_plane._beam = gbt_100m_beam

# gbt_plane_convolved_100 = gbt_plane.convolve_to(ebhis_plane.beam)

# ax = plt.subplot(131)
# hist2d(gbt_plane.value.ravel(), ebhis_plane.value.ravel(), ax=ax)
# plt.plot([0, 15], [0, 15])
# ax2 = plt.subplot(132)
# hist2d(gbt_plane_convolved.value.ravel(), ebhis_plane.value.ravel(), ax=ax2)
# plt.plot([0, 15], [0, 15])
# ax3 = plt.subplot(133)
# hist2d(gbt_plane_convolved_100.value.ravel(), ebhis_plane.value.ravel(), ax=ax3)
# plt.plot([0, 15], [0, 15])


# Best match for GBT is with a 106 m beam, convolved to the 80 m of EBHIS.
# Well, something is wrong here. It has to be that the difference between the
# data is a 80 m deconvolved w/ a 106 m beam. The EBHIS beam size should then
# be slightly smaller?

# Now convolve the Arecibo down to the GBT.
# gbt_90m_beam = Beam(beam_fwhm(90 * u.m))

# arecibo_plane_convolved = arecibo_plane.convolve_to(gbt_90m_beam)


# gbt_100m_beam = Beam(beam_fwhm(100 * u.m))

# arecibo_plane_convolved_100 = arecibo_plane.convolve_to(gbt_100m_beam)

# ax = plt.subplot(131)
# hist2d(arecibo_plane.value.ravel(), gbt_plane.value.ravel(), ax=ax)
# plt.plot([0, 15], [0, 15])
# ax2 = plt.subplot(132)
# hist2d(arecibo_plane_convolved.value.ravel(), gbt_plane.value.ravel(), ax=ax2)
# plt.plot([0, 15], [0, 15])
# ax3 = plt.subplot(133)
# hist2d(arecibo_plane_convolved_100.value.ravel(), gbt_plane.value.ravel(), ax=ax3)
# plt.plot([0, 15], [0, 15])