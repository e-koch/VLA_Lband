
'''
Create a set of thin PV slices
'''

from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy import units as u
import numpy as np

from cube_analysis.disk_pvslices import disk_pvslices

from paths import fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict
from galaxy_params import gal_feath


cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.1kms.fits"))
mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0])


thetas = np.arange(0, 180, 5) * u.deg

pv_width = 40 * u.arcsec

max_rad = 9. * u.kpc

save_name = fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.1kms",
                                        no_check=True)


# Run pv slicing

disk_pvslices(cube, gal_feath, thetas, pv_width, max_rad,
              save_name=save_name, quicklook=False, mom0=mom0,
              save_kwargs=dict(overwrite=True),
              save_regions=True)
