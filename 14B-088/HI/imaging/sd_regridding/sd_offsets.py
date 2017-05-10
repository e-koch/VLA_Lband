
'''
Check the positional offsets between the VLA cube and the SD datasets.
'''

from spectral_cube import SpectralCube
import os
import numpy as np
from astropy import log
from astropy.table import Table

from cube_analysis.register_cubes import cube_registration, spatial_shift_cube

from paths import fourteenB_HI_data_path, data_path

# Load in the 4 cubes and run.

vla_cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))


arecibo_path = os.path.join(data_path, "Arecibo")

# Spectral interpolation, followed by reprojection.
arecibo_name = \
    os.path.join(arecibo_path,
                 "14B-088_items_new/m33_arecibo_14B088_spectralregrid.fits")
arecibo_cube = SpectralCube.read(arecibo_name)

ebhis_path = os.path.join(data_path, "EBHIS")

# Spectral interpolation, followed by reprojection.
ebhis_name = os.path.join(ebhis_path, "14B-088_items/m33_ebhis_14B088_spectralregrid.fits")
ebhis_cube = SpectralCube.read(ebhis_name)


gbt_path = os.path.join(data_path, "GBT")
gbt_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid.fits")
gbt_cube = SpectralCube.read(gbt_name)

gbt_lowres_name = os.path.join(gbt_path, "14B-088_items/m33_gbt_vlsr_Tmb_14B088_spectralregrid.fits")
gbt_lowres_cube = SpectralCube.read(gbt_lowres_name)


num_cores = 6

# Now find the offsets for each
arecibo_offset = cube_registration(arecibo_cube, vla_cube, num_cores=num_cores)
median_arecibo_offset = np.median(arecibo_offset, axis=0)
log.info("Arecibo offset: {0}, {1}".format(*median_arecibo_offset))

ebhis_offset = cube_registration(ebhis_cube, vla_cube, num_cores=num_cores)
median_ebhis_offset = np.median(ebhis_offset, axis=0)
log.info("EBHIS offset: {0}, {1}".format(*median_ebhis_offset))

gbt_offset = cube_registration(gbt_cube, vla_cube, num_cores=num_cores)
median_gbt_offset = np.median(gbt_offset, axis=0)
log.info("GBT offset: {0}, {1}".format(*median_gbt_offset))

gbt_lowres_offset = cube_registration(gbt_lowres_cube, vla_cube,
                                      num_cores=num_cores)
median_gbt_lowres_offset = np.median(gbt_offset, axis=0)
log.info("GBT LR offset: {0}, {1}".format(*median_gbt_lowres_offset))


def save_shift(tablename, pixshifts):
    '''
    Save a csv of the pixel shifts for each channel.
    '''

    t = Table({"y_shift": pixshifts[:, 0], "x_shift": pixshifts[:, 1]})

    t.write(tablename, format='csv')


# Save the pixel shifts
save_shift(os.path.join(arecibo_path, "14B-088_items_new/14B-088_pixel_shifts.csv"))
save_shift(os.path.join(ebhis_path, "14B-088_items/14B-088_pixel_shifts.csv"))
save_shift(os.path.join(gbt_path, "14B-088_items/14B-088_pixel_shifts_highres.csv"))
save_shift(os.path.join(gbt_path, "14B-088_items/14B-088_pixel_shifts.csv"))

# Use the median offsets to shift the cubes

log.info("Shifting and saving Arecibo cube")
arecibo_shift_name = \
    os.path.join(arecibo_path,
                 "14B-088_items_new/m33_arecibo_14B088_spectralregrid_registered.fits")
spatial_shift_cube(arecibo_cube, *median_arecibo_offset, save_shifted=True,
                   save_name=arecibo_shift_name, num_cores=num_cores)

log.info("Shifting and saving EBHIS cube")
ebhis_shift_name = \
    os.path.join(arecibo_path,
                 "14B-088_items/m33_ebhis_14B088_spectralregrid_registered.fits")
spatial_shift_cube(ebhis_cube, *median_ebhis_offset, save_shifted=True,
                   save_name=ebhis_shift_name, num_cores=num_cores)

log.info("Shifting and saving GBT cube")
gbt_shift_name = \
    os.path.join(arecibo_path,
                 "14B-088_items/m33_gbt_vlsr_highres_Tmb_14B088_spectralregrid_registered.fits")
spatial_shift_cube(gbt_cube, *median_gbt_offset, save_shifted=True,
                   save_name=gbt_shift_name, num_cores=num_cores)

log.info("Shifting and saving LR GBT cube")
gbt_lowres_shift_name = \
    os.path.join(arecibo_path,
                 "14B-088_items/m33_gbt_vlsr_Tmb_14B088_spectralregrid_registered.fits")
spatial_shift_cube(gbt_lowres_cube, *median_gbt_lowres_offset,
                   save_shifted=True,
                   save_name=gbt_lowres_shift_name, num_cores=num_cores)
