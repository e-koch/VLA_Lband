
'''
Use the peak velocity diskfit with a radial component to estimate
the radial mass flow in the disk.
'''

from spectral_cube import Projection
from astropy.io import fits
import astropy.units as u
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import seaborn as sb

from constants import hi_mass_conversion, hi_freq, ang_to_phys
from galaxy_params import gal_feath as gal
from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path)
from plotting_styles import onecolumn_figure

from cube_analysis.rotation_curves.radial_inflow import rad_flow

mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0']))

surfdens = (mom0 * (mom0.beam.jtok(hi_freq) * u.beam / u.Jy) *
            hi_mass_conversion).to(u.solMass / u.pc**2) * np.cos(gal.inclination)

tab = Table.read(fourteenB_HI_data_wGBT_path("diskfit_peakvels_noasymm_radial_nowarp_output/rad.out.csv"))

pix_to_phys = ang_to_phys(mom0.header['CDELT2'] * u.deg) / u.pix

mass_flow, mass_flow_err = rad_flow(surfdens, tab, gal, pix_to_phys)

rad_bins = tab['r'] * u.pix * pix_to_phys.to(u.kpc / u.pix)

onecolumn_figure()

plt.errorbar(rad_bins.value, mass_flow.value, yerr=mass_flow_err.value)
plt.grid()
plt.xlabel("R (kpc)")
plt.ylabel(r"Mass Flow Rate (M$_\odot$ yr$^{-1}$)")
plt.axhline(0, linestyle='--', color=sb.color_palette()[1])

plt.tight_layout()

plt.savefig(allfigs_path("rotcurves/radial_mass_flow_rate.png"))
plt.savefig(allfigs_path("rotcurves/radial_mass_flow_rate.pdf"))

plt.close()
