import matplotlib.pyplot as p
import numpy as np
from astropy import units as u
from astropy.stats import mad_std
from spectral_cube import SpectralCube
import seaborn as sb
from corner import hist2d
from os.path import join, exists
import os

from paths import (iram_co21_14B088_data_path, fourteenB_HI_data_path,
                   fourteenB_HI_data_wGBT_path, allfigs_path)
from constants import hi_freq
from galaxy_params import gal
from plotting_styles import default_figure

'''
Regrid the HI data to match the CO.
'''

fig_path = join(allfigs_path(""), "co_vs_hi")
if not exists(fig_path):
    os.mkdir(fig_path)

# Use the CO data smoothed to match the HI
cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.smoothonly.fits"))

hi_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.co21_iram_regrid_full.fits"))
hi_cube_feath = SpectralCube.read(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.co21_iram_regrid_full.fits"))

# Apply masks out to 7 kpc to avoid CO map edges
radii = gal.radius(header=cube.header).to(u.kpc)
max_rad = 6.2 * u.kpc

cube = cube.with_mask(radii < max_rad)
hi_cube = hi_cube.with_mask(radii < max_rad)

co_vals = cube.filled_data[:].flatten()

hi_vals = hi_cube.filled_data[:].flatten()
hi_vals = hi_vals.to(u.K, hi_cube.beam.jtok_equiv(hi_freq))

nan_vals = ~np.logical_or(np.isnan(hi_vals), np.isnan(co_vals))

all_vals = np.vstack([hi_vals.value[nan_vals],
                      co_vals.value[nan_vals]])

# Estimates of the noise level.
# Now the bowls will cause this to be a bit of an overestimation. But we're
# only using this value for the plot. A few far outliers will be down-weighted
# from using the MAD std.
hi_std = mad_std(hi_vals[hi_vals < 0])
co_std = mad_std(co_vals[hi_vals < 0])

default_figure()
sb.set_context("paper", rc={"figure.figsize": np.array([6.4, 5.8])},
               font_scale=1.5)

hist2d(all_vals[0], all_vals[1], bins=50,
       levels=1.0 - np.exp(-0.5 * np.arange(0.5, 3.6, 0.5) ** 2),
       data_kwargs={"alpha": 0.6})
p.hlines(0.0, -20, 110, color='k', linestyle='-')
p.vlines(0.0, -0.2, 1.4, color='k', linestyle='-')

# Noise limits.
p.hlines(3 * co_std.value, -20, 110, color='r', linestyle='--')
p.vlines(3 * hi_std.value, -0.2, 1.4, color='r', linestyle='--')

p.ylim([-0.06, 0.8])

p.ylabel("CO Intensity (K)")
p.xlabel("HI Intensity (K)")

p.savefig(allfigs_path(join(fig_path, "co_vs_hi_brightness_hist2d.png")))
p.savefig(allfigs_path(join(fig_path, "co_vs_hi_brightness_hist2d.pdf")))
p.close()

# Now with the feathered HI
hi_vals_feath = hi_cube_feath.filled_data[:].flatten()
hi_vals_feath = hi_vals_feath.to(u.K, hi_cube.beam.jtok_equiv(hi_freq))

nan_vals = ~np.logical_or(np.isnan(hi_vals_feath), np.isnan(co_vals))

all_vals = np.vstack([hi_vals_feath.value[nan_vals],
                      co_vals.value[nan_vals]])

# Estimates of the noise level.
# Now the bowls will cause this to be a bit of an overestimation. But we're
# only using this value for the plot. A few far outliers will be down-weighted
# from using the MAD std.
hi_std = mad_std(hi_vals_feath[hi_vals_feath < 0])

default_figure()
sb.set_context("paper", rc={"figure.figsize": np.array([6.4, 5.8])},
               font_scale=1.5)

hist2d(all_vals[0], all_vals[1], bins=50,
       levels=1.0 - np.exp(-0.5 * np.arange(0.5, 3.6, 0.5) ** 2),
       data_kwargs={"alpha": 0.6})
p.hlines(0.0, -20, 110, color='k', linestyle='-')
p.vlines(0.0, -0.2, 1.4, color='k', linestyle='-')

# Noise limits.
p.hlines(3 * co_std.value, -20, 110, color='r', linestyle='--')
p.vlines(3 * hi_std.value, -0.2, 1.4, color='r', linestyle='--')

p.ylim([-0.06, 0.8])

p.ylabel("CO Intensity (K)")
p.xlabel("HI Intensity (K)")

p.savefig(allfigs_path(join(fig_path, "co_vs_hi_brightness_hist2d_feath.png")))
p.savefig(allfigs_path(join(fig_path, "co_vs_hi_brightness_hist2d_feath.pdf")))
p.close()

default_figure()
