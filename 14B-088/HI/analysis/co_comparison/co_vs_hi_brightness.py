import matplotlib.pyplot as p
import numpy as np
from astropy import units as u
from astropy.stats import mad_std
from spectral_cube import SpectralCube
import seaborn as sb
# import pandas as pd
from corner import hist2d

from paths import (iram_co21_data_path, fourteenB_HI_data_path,
                   paper1_figures_path)
from constants import regridco_cube_name, hi_freq
from galaxy_params import gal
from plotting_styles import default_figure

'''
Regrid the HI data to match the CO.
'''

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
del cube._header[""]

hi_cube = SpectralCube.read(fourteenB_HI_data_path(regridco_cube_name))

# Apply masks out to 6 kpc to avoid CO map edges
radii = gal.radius(header=cube.header).to(u.kpc)
max_rad = 6 * u.kpc

cube = cube.with_mask(radii < max_rad)
hi_cube = hi_cube.with_mask(radii < max_rad)

co_vals = cube.filled_data[:].flatten()

hi_vals = hi_cube.filled_data[:].flatten()
hi_vals = hi_vals.to(u.K, equivalencies=hi_cube.beam.jtok_equiv(hi_freq))

# df = pd.DataFrame({"CO": co_vals.value, "HI": hi_vals.value})

# sb.jointplot("HI", "CO", data=df, kind='hex')

nan_vals = ~np.logical_or(np.isnan(hi_vals), np.isnan(co_vals))

all_vals = np.vstack([hi_vals.value[nan_vals],
                      co_vals.value[nan_vals]])

# Estimates of the noise level.
# Now the bowls will cause this to be a bit of an overestimation. But we're
# only using this value for the plot. A few far outliers will be down-weighted
# from using the MAD std.
hi_std = mad_std(hi_vals[hi_vals < 0])

default_figure()
sb.set_context("paper", rc={"figure.figsize": np.array([6.4, 5.8])},
               font_scale=1.5)

hist2d(all_vals[0], all_vals[1], bins=50,
       levels=1.0 - np.exp(-0.5 * np.arange(0.5, 3.6, 0.5) ** 2),
       data_kwargs={"alpha": 0.6})
p.hlines(0.0, -20, 110, color='k', linestyle='-')
p.vlines(0.0, -0.2, 1.4, color='k', linestyle='-')

# Noise limits.
p.hlines(3 * 0.02033, -20, 110, color='r', linestyle='--')
p.vlines(3 * hi_std.value, -0.2, 1.4, color='r', linestyle='--')

p.ylabel("CO Intensity (K)")
p.xlabel("HI Intensity (K)")

p.savefig(paper1_figures_path("co_vs_hi_brightness_hist2d.png"))
p.savefig(paper1_figures_path("co_vs_hi_brightness_hist2d.pdf"))
p.close()

default_figure()
