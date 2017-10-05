
from spectral_cube import SpectralCube
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import Angle
import os

from cube_analysis.spectral_stacking import total_profile

from paths import (fourteenB_HI_data_path, fourteenB_HI_file_dict)

from constants import hi_freq

from galaxy_params import gal

hi_cube = SpectralCube.read(fourteenB_HI_file_dict["RotSub_Cube"])
hi_mask = fits.open(fourteenB_HI_file_dict["RotSub_Mask"])[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

hi_cube_cent = SpectralCube.read(fourteenB_HI_file_dict["CentSub_Cube"])
hi_mask_cent = fits.open(fourteenB_HI_file_dict["CentSub_Mask"])[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

hi_cube_peakvel = SpectralCube.read(fourteenB_HI_file_dict["PeakSub_Cube"])
hi_mask_peakvel = fits.open(fourteenB_HI_file_dict["PeakSub_Mask"])[0]
hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask_peakvel.data > 0)

hi_beam = hi_cube.beam

hi_radius = gal.radius(header=hi_cube.header)
hi_pas = gal.position_angles(header=hi_cube.header)

dr = 100 * u.pc
max_radius = (8.0 * u.kpc).to(u.pc)

nbins = np.int(np.floor(max_radius / dr))
inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)

rot_shape = hi_cube.shape[0]
cent_shape = hi_cube_cent.shape[0]
peak_shape = hi_cube_peakvel.shape[0]

total_spectrum_hi_radial = \
    np.zeros((inneredge.size, rot_shape)) * u.K
total_spectrum_hi_radial_n = \
    np.zeros((inneredge.size, rot_shape)) * u.K
total_spectrum_hi_radial_s = \
    np.zeros((inneredge.size, rot_shape)) * u.K
total_spectrum_hi_radial_cent = \
    np.zeros((inneredge.size, cent_shape)) * u.K
total_spectrum_hi_radial_cent_n = \
    np.zeros((inneredge.size, cent_shape)) * u.K
total_spectrum_hi_radial_cent_s = \
    np.zeros((inneredge.size, cent_shape)) * u.K
total_spectrum_hi_radial_peakvel = \
    np.zeros((inneredge.size, peak_shape)) * u.K
total_spectrum_hi_radial_peakvel_n = \
    np.zeros((inneredge.size, peak_shape)) * u.K
total_spectrum_hi_radial_peakvel_s = \
    np.zeros((inneredge.size, peak_shape)) * u.K

# Make the N and S masks
pa_bounds_n = Angle([0.5 * np.pi * u.rad, -0.5 * np.pi * u.rad])
pa_bounds_n = pa_bounds_n.wrap_at(0.5 * np.pi * u.rad)
hi_pa_mask_n = np.logical_and(hi_pas.wrap_at(0.5 * np.pi * u.rad) >= pa_bounds_n[0],
                              hi_pas.wrap_at(0.5 * np.pi * u.rad) < pa_bounds_n[1])

pa_bounds_s = Angle([-0.5 * np.pi * u.rad, 0.5 * np.pi * u.rad])
hi_pa_mask_s = np.logical_and(hi_pas >= pa_bounds_s[0],
                              hi_pas < pa_bounds_s[1])

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {}".format(r0.value, r1))

    hi_rad_mask = np.logical_and(hi_radius >= r0,
                                 hi_radius < r1)

    hi_mask_n = np.logical_and(hi_rad_mask, hi_pa_mask_n)
    hi_mask_s = np.logical_and(hi_rad_mask, hi_pa_mask_s)

    total_spectrum_hi_radial_n[ctr] = \
        total_profile(hi_cube, hi_mask_n).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))
    total_spectrum_hi_radial_s[ctr] = \
        total_profile(hi_cube, hi_mask_s).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_radial[ctr] = \
        total_spectrum_hi_radial_n[ctr] + total_spectrum_hi_radial_s[ctr]

    total_spectrum_hi_radial_cent_n[ctr] = \
        total_profile(hi_cube_cent, hi_mask_n).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))
    total_spectrum_hi_radial_cent_s[ctr] = \
        total_profile(hi_cube_cent, hi_mask_s).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_radial_cent[ctr] = \
        total_spectrum_hi_radial_cent_n[ctr] + \
        total_spectrum_hi_radial_cent_s[ctr]

    total_spectrum_hi_radial_peakvel_n[ctr] = \
        total_profile(hi_cube_peakvel, hi_mask_n).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))
    total_spectrum_hi_radial_peakvel_s[ctr] = \
        total_profile(hi_cube_peakvel, hi_mask_s).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_radial_peakvel[ctr] = \
        total_spectrum_hi_radial_peakvel_n[ctr] + \
        total_spectrum_hi_radial_peakvel_s[ctr]

# Make total versions over the whole cube
total_spectrum_hi = \
    total_profile(hi_cube).to(u.K, hi_beam.jtok_equiv(hi_freq))

total_spectrum_hi_cent = \
    total_profile(hi_cube_cent).to(u.K, hi_beam.jtok_equiv(hi_freq))

total_spectrum_hi_peakvel = \
    total_profile(hi_cube_peakvel).to(u.K, hi_beam.jtok_equiv(hi_freq))


# We'll make mock SpectralCubes from these so it's easy to calculate
# moments
rot_stack_n = SpectralCube(data=total_spectrum_hi_radial_n.T.reshape((rot_shape, inneredge.size, 1)),
                           wcs=hi_cube.wcs)
rot_stack_s = SpectralCube(data=total_spectrum_hi_radial_s.T.reshape((rot_shape, inneredge.size, 1)),
                           wcs=hi_cube.wcs)
rot_stack = SpectralCube(data=total_spectrum_hi_radial.T.reshape((rot_shape, inneredge.size, 1)),
                         wcs=hi_cube.wcs)

cent_stack_n = SpectralCube(data=total_spectrum_hi_radial_cent_n.T.reshape((cent_shape, inneredge.size, 1)),
                            wcs=hi_cube.wcs)
cent_stack_s = SpectralCube(data=total_spectrum_hi_radial_cent_s.T.reshape((cent_shape, inneredge.size, 1)),
                            wcs=hi_cube.wcs)
cent_stack = SpectralCube(data=total_spectrum_hi_radial_cent.T.reshape((cent_shape, inneredge.size, 1)),
                          wcs=hi_cube.wcs)

peakvel_stack_n = SpectralCube(data=total_spectrum_hi_radial_peakvel_n.T.reshape((peak_shape, inneredge.size, 1)),
                               wcs=hi_cube.wcs)
peakvel_stack_s = SpectralCube(data=total_spectrum_hi_radial_peakvel_s.T.reshape((peak_shape, inneredge.size, 1)),
                               wcs=hi_cube.wcs)
peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((peak_shape, inneredge.size, 1)),
                             wcs=hi_cube.wcs)

# Now save all of these for future use.
stacked_folder = fourteenB_HI_data_path("stacked_spectra", no_check=True)
if not os.path.exists(stacked_folder):
    os.mkdir(stacked_folder)

wstring = "{0}{1}".format(int(dr.value), dr.unit)
rot_stack_n.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_north_{}.fits".format(wstring),
                                         no_check=True), overwrite=True)
rot_stack_s.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_south_{}.fits".format(wstring),
                                         no_check=True), overwrite=True)
rot_stack.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_radial_{}.fits".format(wstring),
                                       no_check=True), overwrite=True)

cent_stack_n.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_north_{}.fits".format(wstring),
                                          no_check=True), overwrite=True)
cent_stack_s.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_south_{}.fits".format(wstring),
                                          no_check=True), overwrite=True)
cent_stack.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_radial_{}.fits".format(wstring),
                                        no_check=True), overwrite=True)

peakvel_stack_n.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring),
                                             no_check=True), overwrite=True)
peakvel_stack_s.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring),
                                             no_check=True), overwrite=True)
peakvel_stack.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring),
                                           no_check=True), overwrite=True)

total_spectrum_hi.hdu.writeto(fourteenB_HI_data_path("stacked_spectra/rotation_stacked.fits",
                                                     no_check=True),
                              overwrite=True)
total_spectrum_hi.hdu.writeto(fourteenB_HI_data_path("stacked_spectra/rotation_stacked.fits",
                                                     no_check=True),
                              overwrite=True)
total_spectrum_hi.hdu.writeto(fourteenB_HI_data_path("stacked_spectra/rotation_stacked.fits",
                                                     no_check=True),
                              overwrite=True)
