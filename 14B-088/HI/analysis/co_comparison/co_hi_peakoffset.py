
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as p
from scipy import ndimage as nd
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.utils.console import ProgressBar

from rotation_curves.vrot_fit import return_smooth_model
from paths import (fourteenB_HI_data_path, iram_co21_data_path,
                   paper1_figures_path)
from galaxy_params import gal

try:
    from corner import hist2d
except ImportError:
    raise ImportError("You need to install the corner package!")

'''
Is bright CO better correlated w/ the HI velocity?

Compares Tpeak_CO vs |V_peak,CO - V_rot,HI|
'''

table_name = fourteenB_HI_data_path("diskfit_noasymm_noradial_nowarp_output/rad.out.csv")
tab = Table.read(table_name)

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
del cube._header[""]

# Galactic radii
radii = gal.radius(header=cube.header)

hi_model = return_smooth_model(tab, cube.header, gal) * u.m / u.s

# Mask everything at the edges of the map.
r_max = 6. * u.kpc
hi_model[radii > r_max] = np.NaN * u.m / u.s
cube = cube.with_mask(radii < r_max)

# Use avg noise reported in Druard+14
co_avg_noise = 20.33 * u.mK

# Only keep peaks above 3 sigma.
min_snr = 5.

Tpeak = np.zeros(cube.shape[1:]) * u.K
vel_Tpeak = np.zeros(cube.shape[1:]) * u.m / u.s

good_posns = np.where(cube.mask.include().sum(0) > 0)

for y, x in ProgressBar(zip(*good_posns)):
    specind = cube[:, y, x].argmax()
    Tpeak[y, x] = cube[:, y, x][specind]

    if Tpeak[y, x] < min_snr * co_avg_noise:
        Tpeak[y, x] = np.NaN * u.K
        vel_Tpeak[y, x] = np.NaN * u.m / u.s
    else:
        vel_Tpeak[y, x] = cube.spectral_axis[specind]

Tpeak[Tpeak == 0.0 * u.K] = np.NaN * u.K
vel_Tpeak[vel_Tpeak == 0.0 * u.m / u.s] = np.NaN * u.m / u.s

vel_diff = np.abs(hi_model - vel_Tpeak).to(u.km / u.s)

# Clean-up the points by requiring they be part of a region larger than the
# beam
Tpeak_mask = np.isfinite(Tpeak)
pixscale = proj_plane_pixel_scales(cube.wcs)[0]
beam_kern = cube.beam.as_tophat_kernel(pixscale)

# Open and close
Tpeak_mask = nd.binary_opening(Tpeak_mask, beam_kern)
Tpeak_mask = nd.binary_closing(Tpeak_mask, beam_kern)

vel_diff_values = vel_diff.value[Tpeak_mask & np.isfinite(Tpeak)]
Tpeak_values = Tpeak.value[Tpeak_mask & np.isfinite(Tpeak)]

hist2d(vel_diff_values, Tpeak_values, bins=16, data_kwargs={"alpha": 0.6},
       range=[(0.0, 30.0),
              (0.0, np.max(Tpeak_values))])
p.hlines(co_avg_noise.to(u.K).value * min_snr, 0.0, 30.0, color='r',
         linestyle='--')
p.ylabel(r"T$_\mathrm{peak, CO}$ (K)")
p.xlabel(r"|V$_\mathrm{peak, CO}$ - V$_\mathrm{rot, HI}$| (km/s)")
p.grid()

p.savefig(paper1_figures_path("co21_Tpeak_velocity_offset.pdf"))
p.savefig(paper1_figures_path("co21_Tpeak_velocity_offset.png"))
p.close()
