
'''
Exploring the properties of the skewness and kurtosis arrays.
'''

from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import numpy as np
import matplotlib.pyplot as p
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle

from analysis.paths import (fourteenB_HI_data_path, paper1_figures_path,
                            iram_co21_data_path, data_path)
from constants import (hi_freq, cube_name, moment0_name, lwidth_name,
                       skew_name, kurt_name, mask_name, moment1_name)
from plotting_styles import (twocolumn_figure, onecolumn_figure,
                             default_figure, twocolumn_twopanel_figure)

from HI_veldisp_profile import radial_profile
from galaxy_params import gal


mom0_hdu = fits.open(fourteenB_HI_data_path(moment0_name))[0]
mom0 = Projection.from_hdu(mom0_hdu)

mom1_hdu = fits.open(fourteenB_HI_data_path(moment1_name))[0]
mom1 = Projection.from_hdu(mom1_hdu)

lwidth_hdu = fits.open(fourteenB_HI_data_path(lwidth_name))[0]
lwidth = Projection.from_hdu(lwidth_hdu)

skew_hdu = fits.open(fourteenB_HI_data_path(skew_name))[0]
# Remove a few bad outliers from skew
skew_hdu.data[skew_hdu.data > 10] = np.NaN
skew_hdu.data[skew_hdu.data < -10] = np.NaN
skew = Projection.from_hdu(skew_hdu)

kurt_hdu = fits.open(fourteenB_HI_data_path(kurt_name))[0]
kurt_hdu.data[kurt_hdu.data > 10] = np.NaN
kurt_hdu.data[kurt_hdu.data < -10] = np.NaN
kurt = Projection.from_hdu(kurt_hdu)


dr = 100 * u.pc

radii_skew, sdprof_skew, sdprof_sigma_skew = \
    radial_profile(gal, skew, dr=dr, max_rad=8 * u.kpc)

rs_skew_n, sd_skew_n, sd_skew_sigma_n = \
    radial_profile(gal, skew, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
rs_skew_s, sd_skew_s, sd_skew_sigma_s = \
    radial_profile(gal, skew, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))

radii_kurt, sdprof_kurt, sdprof_sigma_kurt = \
    radial_profile(gal, kurt, dr=dr, max_rad=8 * u.kpc)

rs_kurt_n, sd_kurt_n, sd_kurt_sigma_n = \
    radial_profile(gal, kurt, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([0.5 * np.pi * u.rad,
                                    -0.5 * np.pi * u.rad]))
rs_kurt_s, sd_kurt_s, sd_kurt_sigma_s = \
    radial_profile(gal, kurt, max_rad=8 * u.kpc, dr=dr,
                   pa_bounds=Angle([-0.5 * np.pi * u.rad,
                                    0.5 * np.pi * u.rad]))

twocolumn_twopanel_figure()

p.subplot(121)
# p.errorbar(radii.value, sdprof.value, yerr=sdprof_sigma.value,
#            drawstyle='steps-mid', label="Total", color='b', zorder=-1)
p.plot(radii_skew.value, sdprof_skew.value,
       drawstyle='steps-mid', label="Total", color='b', zorder=-1)
p.plot(rs_skew_n.value, sd_skew_n.value, 'r-.', drawstyle='steps-mid',
       label='North', alpha=0.5)
p.plot(rs_skew_s.value, sd_skew_s.value, 'g--', drawstyle='steps-mid',
       label='South', alpha=0.5)
# p.errorbar(rs_skew_n.value, sd_skew_n.value, yerr=sd_skew_sigma_n.value,
#            drawstyle='steps-mid', label="North", alpha=0.5)
# p.errorbar(rs_s.value, sd_s.value, yerr=sd_sigma_s.value,
#            drawstyle='steps-mid', label="South", alpha=0.5)
p.legend(frameon=True, loc='lower right')
p.grid()
p.xlabel("Radius (kpc)")
p.ylabel("Skewness")

p.subplot(122)
# p.errorbar(radii_kurt.value, sdprof_kurt.value, yerr=sdprof_sigma_kurt.value,
#            drawstyle='steps-mid', label="Total", color='b', zorder=-1)
p.plot(radii_kurt.value, sdprof_kurt.value,
       drawstyle='steps-mid', label="Total", color='b', zorder=-1)
p.plot(rs_kurt_n.value, sd_kurt_n.value, 'r-.', drawstyle='steps-mid',
       label='North', alpha=0.5)
p.plot(rs_kurt_s.value, sd_kurt_s.value, 'g--', drawstyle='steps-mid',
       label='South', alpha=0.5)
# p.errorbar(rs_kurt_n.value, sd_kurt_n.value, yerr=sd_kurt_sigma_n.value,
#            drawstyle='steps-mid', label="North", alpha=0.5)
# p.errorbar(rs_s.value, sd_s.value, yerr=sd_sigma_s.value,
#            drawstyle='steps-mid', label="South", alpha=0.5)
p.grid()
p.xlabel("Radius (kpc)")
p.ylabel("Kurtosis")

p.tight_layout()

p.savefig(paper1_figures_path("hi_skew_kurt_profile_n_s.png"))
p.savefig(paper1_figures_path("hi_skew_kurt_profile_n_s.pdf"))
p.close()

default_figure()
