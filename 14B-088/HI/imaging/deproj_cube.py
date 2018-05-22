
'''
Create a deprojected cube in M33's frame
'''

from spectral_cube import SpectralCube
from astropy.io import fits
import numpy as np
import astropy.units as u
from radio_beam import Beam

from cube_analysis.cube_deproject import deproject_cube

from paths import (fourteenB_wGBT_HI_file_dict, allfigs_path,
                   fourteenB_HI_data_wGBT_path, data_path)
from galaxy_params import gal_feath as gal

cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.com_beam.fits"))
deproject_cube(cube, gal, num_cores=6, chunk=100,
               save_name=fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.deproject.fits",
                                                     no_check=True))

hdu = fits.open(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.deproject.fits"),
                mode='update')

# Update the beam in the header
hdr = hdu[0].header
hdr.update(Beam(major=cube.beam.major / np.cos(gal.inclination),
                minor=cube.beam.major,
                pa=gal.position_angle + 90 * u.deg).to_header_keywords())
hdu[0].header = hdr

hdu.flush()
hdu.close()

# Do the same for the peak-velocity corrected cube
cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peakvels_corrected.fits"))
# Limit the velocity range to where the actual emission is. Ranges set by-eye.
deproject_cube(cube[595:1372], gal, num_cores=6, chunk=40,
               save_name=fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peakvels_corrected.deproject.fits",
                                                     no_check=True))

hdu = fits.open(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peakvels_corrected.deproject.fits"),
                mode='update')

# Update the beam in the header
hdr = hdu[0].header
hdr.update(Beam(major=cube.beam.major / np.cos(gal.inclination),
                minor=cube.beam.major,
                pa=gal.position_angle + 90 * u.deg).to_header_keywords())
hdu[0].header = hdr

hdu.flush()
hdu.close()
