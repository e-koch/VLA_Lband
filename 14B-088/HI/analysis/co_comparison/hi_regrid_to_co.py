
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import largest_beam

from paths import iram_co21_data_path, fourteenB_HI_data_path
from constants import cube_name, regridco_cube_name

'''
Regrid the HI data to match the CO.
'''

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
del cube._header[""]

hi_cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

# Find the largest HI beam
hi_beam = largest_beam(hi_cube.beams)

hi_conv_cube = hi_cube.convolve_to(hi_beam)

hi_conv_cube = hi_conv_cube.spectral_interpolate(cube.spectral_axis)

hi_conv_cube = hi_conv_cube.reproject(cube.header)

# For some reason, the velocity axis in the data needs to be flipped.
# Look into where this is happening in spectral_cube
hi_hdu = hi_conv_cube.hdu.copy()
hi_hdu.data = hi_hdu.data[::-1]

new_hi_conv_cube = SpectralCube.read(hi_hdu)

new_hi_conv_cube.write(fourteenB_HI_data_path(regridco_cube_name, no_check=True))
