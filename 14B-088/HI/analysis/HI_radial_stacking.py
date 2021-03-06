
from spectral_cube import SpectralCube
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import Angle
import os

from cube_analysis.spectral_stacking import total_profile, radial_stacking

from paths import (fourteenB_HI_data_path, fourteenB_HI_file_dict)

from constants import hi_freq

from galaxy_params import gal

hi_cube = SpectralCube.read(fourteenB_HI_file_dict["RotSub_Cube"],
                            memmap=True)

hi_beam = hi_cube.beam
jybm_to_K = hi_beam.jtok_equiv(hi_freq)

hi_radius = gal.radius(header=hi_cube.header)
hi_pas = gal.position_angles(header=hi_cube.header)

dr = 100 * u.pc
max_radius = (8.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)

nbins = np.int(np.floor(max_radius / dr))
inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)
# Make the N and S masks
pa_bounds_n = Angle([0.5 * np.pi * u.rad, -0.5 * np.pi * u.rad])

pa_bounds_s = Angle([-0.5 * np.pi * u.rad, 0.5 * np.pi * u.rad])

hi_mask = fits.open(fourteenB_HI_file_dict["RotSub_Mask"])[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

del hi_mask

rot_shape = hi_cube.shape[0]

total_spectrum_hi_radial = radial_stacking(gal, hi_cube, dr=dr,
                                           max_radius=max_radius,
                                           pa_bounds=None,
                                           verbose=True,
                                           how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_radial_n = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_radial_s = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

# Total over all spectra
total_spectrum_hi = \
    total_profile(hi_cube).to(u.K, jybm_to_K)

# Now save all of these for future use.
stacked_folder = fourteenB_HI_data_path("stacked_spectra", no_check=True)
if not os.path.exists(stacked_folder):
    os.mkdir(stacked_folder)

# Radial stacks
rot_stack_n = SpectralCube(data=total_spectrum_hi_radial_n.T.reshape((rot_shape, inneredge.size, 1)),
                           wcs=hi_cube.wcs)
rot_stack_s = SpectralCube(data=total_spectrum_hi_radial_s.T.reshape((rot_shape, inneredge.size, 1)),
                           wcs=hi_cube.wcs)
rot_stack = SpectralCube(data=total_spectrum_hi_radial.T.reshape((rot_shape, inneredge.size, 1)),
                         wcs=hi_cube.wcs)

rot_stack_n.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_north_{}.fits".format(wstring),
                                         no_check=True), overwrite=True)
rot_stack_s.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_south_{}.fits".format(wstring),
                                         no_check=True), overwrite=True)
rot_stack.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_{}.fits".format(wstring),
                                       no_check=True), overwrite=True)

total_spectrum_hi.hdu.writeto(fourteenB_HI_data_path("stacked_spectra/rotation_stacked.fits",
                                                     no_check=True),
                              overwrite=True)

# Now delete that cube to save memory before loading the next
del hi_cube

# Centroid stack
hi_cube_cent = SpectralCube.read(fourteenB_HI_file_dict["CentSub_Cube"],
                                 memmap=True)
hi_mask_cent = fits.open(fourteenB_HI_file_dict["CentSub_Mask"])[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

del hi_mask_cent

cent_shape = hi_cube_cent.shape[0]

total_spectrum_hi_radial_cent = \
    radial_stacking(gal, hi_cube_cent, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_radial_cent_n = \
    radial_stacking(gal, hi_cube_cent, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_radial_cent_s = \
    radial_stacking(gal, hi_cube_cent, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_cent = \
    total_profile(hi_cube_cent).to(u.K, jybm_to_K)

cent_stack_n = SpectralCube(data=total_spectrum_hi_radial_cent_n.T.reshape((cent_shape, inneredge.size, 1)),
                            wcs=hi_cube_cent.wcs)
cent_stack_s = SpectralCube(data=total_spectrum_hi_radial_cent_s.T.reshape((cent_shape, inneredge.size, 1)),
                            wcs=hi_cube_cent.wcs)
cent_stack = SpectralCube(data=total_spectrum_hi_radial_cent.T.reshape((cent_shape, inneredge.size, 1)),
                          wcs=hi_cube_cent.wcs)

cent_stack_n.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_north_{}.fits".format(wstring),
                                          no_check=True), overwrite=True)
cent_stack_s.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_south_{}.fits".format(wstring),
                                          no_check=True), overwrite=True)
cent_stack.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_{}.fits".format(wstring),
                                        no_check=True), overwrite=True)


total_spectrum_hi_cent.hdu.writeto(fourteenB_HI_data_path("stacked_spectra/centroid_stacked.fits",
                                                               no_check=True),
                                   overwrite=True)

del hi_cube_cent


hi_cube_peakvel = SpectralCube.read(fourteenB_HI_file_dict["PeakSub_Cube"])
hi_mask_peakvel = fits.open(fourteenB_HI_file_dict["PeakSub_Mask"])[0]
hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask_peakvel.data > 0)

del hi_mask_peakvel

peak_shape = hi_cube_peakvel.shape[0]


total_spectrum_hi_radial_peakvel = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_radial_peakvel_n = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_radial_peakvel_s = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='cube')[1].to(u.K, jybm_to_K)

total_spectrum_hi_peakvel = \
    total_profile(hi_cube_peakvel).to(u.K, jybm_to_K)


peakvel_stack_n = SpectralCube(data=total_spectrum_hi_radial_peakvel_n.T.reshape((peak_shape, inneredge.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack_s = SpectralCube(data=total_spectrum_hi_radial_peakvel_s.T.reshape((peak_shape, inneredge.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((peak_shape, inneredge.size, 1)),
                             wcs=hi_cube_peakvel.wcs)

peakvel_stack_n.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring),
                                                  no_check=True), overwrite=True)
peakvel_stack_s.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring),
                                                  no_check=True), overwrite=True)
peakvel_stack.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring),
                                                no_check=True), overwrite=True)

total_spectrum_hi_peakvel.hdu.writeto(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked.fits",
                                                                  no_check=True),
                                      overwrite=True)

del hi_cube_peakvel
