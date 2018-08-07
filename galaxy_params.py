
'''
Use parameters from Diskfit in the Galaxy class
'''

from galaxies import Galaxy
from astropy.table import Table
import astropy.units as u

from cube_analysis.rotation_curves import update_galaxy_params

from paths import fourteenB_HI_data_path, fourteenB_HI_data_wGBT_path

# The models from the peak velocity aren't as biased, based on comparing
# the VLA and VLA+GBT velocity curves. Using these as the defaults

folder_name = "diskfit_peakvels_noasymm_noradial_nowarp_output"

param_name = \
    fourteenB_HI_data_path("{}/rad.out.params.csv".format(folder_name))

param_table = Table.read(param_name)

gal = Galaxy("M33")

update_galaxy_params(gal, param_table)

# Load in the model from the feathered data as well.
folder_name = "diskfit_peakvels_noasymm_noradial_nowarp_output"

param_name = \
    fourteenB_HI_data_wGBT_path("{}/rad.out.params.csv".format(folder_name))

param_table = Table.read(param_name)

gal_feath = Galaxy("M33")

update_galaxy_params(gal_feath, param_table)

# Force 840 kpc for the distance

gal.distance = 840 * u.kpc
gal_feath.distance = 840 * u.kpc
