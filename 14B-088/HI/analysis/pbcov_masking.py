
'''
Cut out noisy regions by imposing a mask of the primary beam coverage.
'''

from astropy.io import fits
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import beams_to_bintable
from astropy.utils.console import ProgressBar
import os

from analysis.paths import fourteenB_HI_data_path
from analsis.constants import pb_lim, cube_name

# execfile(os.path.expanduser("~/Dropbox/code_development/ewky_scripts/write_huge_fits.py"))

pbcov = fits.open(fourteenB_HI_data_path("M33_14B-088_pbcov.fits"))[0]

cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.fits"))

# Apply the mask, using a cut-off of 0.3. This retains all of the regions with
# emission.
masked_cube = cube.with_mask(pbcov.data > pb_lim)

masked_cube = masked_cube.minimal_subcube()

new_fitsname = fourteenB_HI_data_path(cube_name, no_check=True)

masked_cube.write(new_fitsname)

# create_huge_fits(new_fitsname, cube.header)

# save_hdu = fits.open(new_fitsname, mode='update')

# Save per channel
# for chan in ProgressBar(cube.shape[0]):
#     save_hdu[0].data[chan] = cube[chan].value

#     if chan % 50 == 0:
#         save_hdu.flush()

# Save the beam table!
# save_hdu.append(beams_to_bintable(cube.beams))

# save_hdu.flush()
# save_hdu.close()
