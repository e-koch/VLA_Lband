
'''
Make a version of the CO(2-1) data at the resolution of the 14B-088 HI data
'''

from spectral_cube import SpectralCube, Projection
from spectral_cube.cube_utils import largest_beam
from astropy.io import fits


from paths import fourteenB_HI_file_dict, iram_co21_data_path


vla_cube = SpectralCube.read(fourteenB_HI_file_dict["Cube"])

co21_cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))

# Smooth to the largest beam (they differ by tiny fractions anyway)

large_beam = largest_beam(vla_cube.beams)
smoothed_cube = co21_cube.convolve_to(large_beam)
smoothed_cube.write(iram_co21_data_path("m33.co21_iram.HI_14B_smoothed.fits", no_check=True))

# Now smooth the noise map
co21_noise = Projection.from_hdu(fits.open(iram_co21_data_path("m33.rms.fits"))[0])
co21_noise_smoothed = co21_noise.convolve_to(large_beam)
co21_noise_smoothed.write(iram_co21_data_path("m33.rms.HI_14B_smoothed.fits", no_check=True))

