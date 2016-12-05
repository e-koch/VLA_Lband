
'''
Use parameters from Diskfit in the Galaxy class
'''

from astropy import units as u
from galaxies import Galaxy
from astropy.table import Table

from paths import fourteenB_HI_data_path


def update_galaxy_params(gal, param_table):
    '''
    Use the fit values from fit rather than the hard-coded values in galaxies.
    '''

    from astropy.coordinates import Angle, SkyCoord

    gal.inclination = Angle(param_table["inc"] * u.deg)[0]
    gal.position_angle = Angle(param_table["PA"] * u.deg)[0]
    gal.vsys = (param_table["Vsys"] * u.km / u.s)[0]

    # The positions in the table are in pixels, so convert to the sky using
    # the spatial WCS info.
    ra_cent, dec_cent = param_table["RAcent"], param_table["Deccent"]

    gal.center_position = SkyCoord(ra_cent, dec_cent, unit=(u.deg, u.deg),
                                   frame='fk5')


folder_name = "diskfit_noasymm_noradial_nowarp_output"

param_name = \
    fourteenB_HI_data_path("{}/rad.out.params.csv".format(folder_name))

param_table = Table.read(param_name)

gal = Galaxy("M33")

update_galaxy_params(gal, param_table)
