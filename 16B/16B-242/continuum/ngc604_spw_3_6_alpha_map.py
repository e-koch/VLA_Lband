
'''
Calculate spectral slope of bright sources in NGC 604.
'''

import numpy as np
from astropy.io import fits
from astropy import units as u
from spectral_cube import Projection
import matplotlib.pyplot as plt

spatial_slice = (slice(3800, 4200),) * 2

# spw3 = Projection.from_hdu(fits.open("NGC604_16B-242_continuum_spw_3_A_B_C.image.fits"))
# spw6 = Projection.from_hdu(fits.open("NGC604_16B-242_continuum_spw_6_A_B_C.image.fits"))

spw3 = Projection.from_hdu(fits.open("NGC604_16B-242_continuum_spw_3_A.image.fits"))
spw6 = Projection.from_hdu(fits.open("NGC604_16B-242_continuum_spw_6_A.image.fits"))

alp_map = Projection.from_hdu(fits.open("NGC604_16B-242_continuum_spw_3_6_A_mtmfs.alpha.fits"))
alperr_map = Projection.from_hdu(fits.open("NGC604_16B-242_continuum_spw_3_6_A_mtmfs.alpha.error.fits"))
mfs_map = Projection.from_hdu(fits.open("NGC604_16B-242_continuum_spw_3_6_A_mtmfs.image.tt0.fits"))

spw3 = spw3[spatial_slice]
spw6 = spw6[spatial_slice]

alp_map = alp_map[spatial_slice]
alperr_map = alperr_map[spatial_slice]
mfs_map = mfs_map[spatial_slice]

# Convolve spw 6 to the 3 beam
spw6 = spw6.convolve_to(spw3.beam)

spw3_freq = 1.472 * u.GHz + 64 * u.MHz
spw6_freq = 1.756 * u.GHz + 64 * u.MHz


alpha_map = ((np.log10(spw6.value) - np.log10(spw3.value)) / (np.log10(spw6_freq.value) - np.log10(spw3_freq.value)))

int_cutoff = 0.0001
# int_cutoff = 0.0002

alpha_map[spw3.value < int_cutoff] = np.NaN


# plt.subplot(projection=spw3.wcs)
# # plt.imshow(np.arcsinh(spw3.value / np.nanpercentile(spw3.value, 95)), origin='lower', cmap='binary', alpha=0.5, zorder=-1, vmin=0.)
# plt.imshow(spw3.value, origin='lower', cmap='binary_r', alpha=0.7, zorder=-1, vmax=0.0005, vmin=0.)
# # plt.imshow(spw6.value, origin='lower', cmap='binary_r', alpha=0.5, zorder=-1, vmax=0.0005, vmin=0.)
# plt.imshow(alpha_map, origin='lower', cmap='viridis', alpha=0.6, vmax=-1.5, vmin=-3.5)
# plt.colorbar()
# cs = plt.contour(alpha_map, cmap='viridis', levels=np.round(np.linspace(-4, 0, 20), 2))
# cs = plt.contour(alpha_map, cmap='viridis', levels=[-4., -3., -2., -1])
# plt.colorbar(cs)

# plt.imshow(alpha_map, origin='lower', cmap='inferno', vmin=-4.0, vmax=-0.5)
# plt.colorbar()
# plt.contour(spw3.value, cmap='viridis')

plt.subplot(projection=mfs_map.wcs)
plt.imshow(mfs_map.value, origin='lower', cmap='binary_r', alpha=0.7, zorder=-1, vmax=0.0005, vmin=0.)
# plt.imshow(spw6.value, origin='lower', cmap='binary_r', alpha=0.5, zorder=-1, vmax=0.0005, vmin=0.)

# alpha_mask = (np.abs(alp_map) / alperr_map) > 3.
alpha_mask = alperr_map < 0.5

alp_map[~alpha_mask] = np.NaN

plt.imshow(alpha_map, origin='lower', cmap='viridis', alpha=0.6, vmax=-1.5, vmin=-3.5)
plt.colorbar()