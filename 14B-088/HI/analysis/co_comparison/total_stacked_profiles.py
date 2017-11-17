
import astropy.units as u
from spectral_cube import SpectralCube, OneDSpectrum
import numpy as np
from astropy.io import fits
from os.path import join as osjoin
import os
from astropy.coordinates import Angle

from cube_analysis.spectral_stacking import radial_stacking

from paths import (fourteenB_wGBT_HI_file_dict,
                   iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path)
from constants import hi_freq
from galaxy_params import gal_feath as gal

'''
Create profiles of HI and CO after subtracting velocity surfaces
'''

co_stackpath = lambda x: osjoin(iram_co21_14B088_data_path("", no_check=True), "stacked_spectra", x)
if not os.path.exists(co_stackpath("")):
    os.mkdir(co_stackpath(""))

dr = 500 * u.pc
max_radius = (7.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)
maxrad_string = "{0}{1}".format(int(max_radius.value), max_radius.unit)

pa_bounds_n = Angle([0.5 * np.pi * u.rad, -0.5 * np.pi * u.rad])

pa_bounds_s = Angle([-0.5 * np.pi * u.rad, 0.5 * np.pi * u.rad])

# CO stacking first
co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_feather.rotation_corrected.fits"))

bin_centers, total_spectrum_co_radial, num_pixels = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='cube')

bin_centers, total_spectrum_co_radial_n, num_pixels_n = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='cube')

bin_centers, total_spectrum_co_radial_s, num_pixels_s = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='cube')

spec_shape = co_cube.shape[0]

rot_stack = SpectralCube(data=total_spectrum_co_radial.T.reshape((spec_shape, bin_centers.size, 1)),
                         wcs=co_cube.wcs)
rot_stack.write(co_stackpath("rotation_stacked_radial_{}.fits".format(wstring)),
                overwrite=True)
rot_stack_n = SpectralCube(data=total_spectrum_co_radial_n.T.reshape((spec_shape, bin_centers.size, 1)),
                           wcs=co_cube.wcs)
rot_stack_n.write(co_stackpath("rotation_stacked_radial_north_{}.fits".format(wstring)),
                  overwrite=True)
rot_stack_s = SpectralCube(data=total_spectrum_co_radial_s.T.reshape((spec_shape, bin_centers.size, 1)),
                           wcs=co_cube.wcs)
rot_stack_s.write(co_stackpath("rotation_stacked_radial_south_{}.fits".format(wstring)),
                  overwrite=True)

# Separately save the number of pixels in each bin
np.save(co_stackpath("radial_stacking_pixelsinbin_{}.npy").format(wstring), num_pixels)
np.save(co_stackpath("radial_stacking_pixelsinbin_north_{}.npy").format(wstring), num_pixels_n)
np.save(co_stackpath("radial_stacking_pixelsinbin_south_{}.npy").format(wstring), num_pixels_s)


# Save the total profiles over the inner 7 kpc

total_spectrum_co = total_spectrum_co_radial.sum(0)

oned_wcs = co_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co.value, unit=total_spectrum_co.unit,
             wcs=oned_wcs).write(co_stackpath("rotation_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del co_cube

co_cube_cent = \
    SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_feather.centroid_corrected.fits"))

total_spectrum_co_radial_cent = \
    radial_stacking(gal, co_cube_cent, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='cube')[1]

total_spectrum_co_radial_cent_n = \
    radial_stacking(gal, co_cube_cent, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='cube')[1]

total_spectrum_co_radial_cent_s = \
    radial_stacking(gal, co_cube_cent, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='cube')[1]

spec_shape = co_cube_cent.shape[0]


cent_stack = SpectralCube(data=total_spectrum_co_radial_cent.T.reshape((spec_shape, bin_centers.size, 1)),
                          wcs=co_cube_cent.wcs)
cent_stack.write(co_stackpath("centroid_stacked_radial_{}.fits".format(wstring)),
                 overwrite=True)

cent_stack_n = SpectralCube(data=total_spectrum_co_radial_cent_n.T.reshape((spec_shape, bin_centers.size, 1)),
                            wcs=co_cube_cent.wcs)
cent_stack_n.write(co_stackpath("centroid_stacked_radial_north_{}.fits".format(wstring)),
                   overwrite=True)
cent_stack_s = SpectralCube(data=total_spectrum_co_radial_cent_s.T.reshape((spec_shape, bin_centers.size, 1)),
                            wcs=co_cube_cent.wcs)
cent_stack_s.write(co_stackpath("centroid_stacked_radial_south_{}.fits".format(wstring)),
                   overwrite=True)

total_spectrum_co_cent = total_spectrum_co_radial_cent.sum(0)
oned_wcs = co_cube_cent[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co_cent.value, unit=total_spectrum_co_cent.unit,
             wcs=oned_wcs).write(co_stackpath("centroid_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del co_cube_cent

co_cube_peakvel = \
    SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_feather.peakvels_corrected.fits"))

total_spectrum_co_radial_peakvel = \
    radial_stacking(gal, co_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='cube')[1]

total_spectrum_co_radial_peakvel_n = \
    radial_stacking(gal, co_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='cube')[1]

total_spectrum_co_radial_peakvel_s = \
    radial_stacking(gal, co_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='cube')[1]

spec_shape = co_cube_peakvel.shape[0]

peakvel_stack = SpectralCube(data=total_spectrum_co_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                          wcs=co_cube_peakvel.wcs)
peakvel_stack.write(co_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)),
                    overwrite=True)

peakvel_stack_n = SpectralCube(data=total_spectrum_co_radial_peakvel_n.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=co_cube_peakvel.wcs)
peakvel_stack_n.write(co_stackpath("peakvel_stacked_radial_north_{}.fits".format(wstring)),
                      overwrite=True)
peakvel_stack_s = SpectralCube(data=total_spectrum_co_radial_peakvel_s.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=co_cube_peakvel.wcs)
peakvel_stack_s.write(co_stackpath("peakvel_stacked_radial_south_{}.fits".format(wstring)),
                      overwrite=True)

total_spectrum_co_peakvel = total_spectrum_co_radial_peakvel.sum(0)

# Save each of these
oned_wcs = co_cube_peakvel[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co_peakvel.value,
             unit=total_spectrum_co_peakvel.unit,
             wcs=oned_wcs).write(co_stackpath("peakvel_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del co_cube_peakvel

# Now do HI stacking with 500 pc bins

hi_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["RotSub_Cube"])
hi_mask = fits.open(fourteenB_wGBT_HI_file_dict["RotSub_Mask"])[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

del hi_mask

hi_beam = hi_cube.beam
jybm_to_K = hi_beam.jtok_equiv(hi_freq)

total_spectrum_hi_radial, num_pixels = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=True,
                    how='cube')[1:]

total_spectrum_hi_radial_n, num_pixels_n = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_n,
                    verbose=True,
                    how='cube')[1:]

total_spectrum_hi_radial_s, num_pixels_s = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=pa_bounds_s,
                    verbose=True,
                    how='cube')[1:]

total_spectrum_hi_radial = total_spectrum_hi_radial.to(u.K, jybm_to_K)
total_spectrum_hi_radial_n = total_spectrum_hi_radial_n.to(u.K, jybm_to_K)
total_spectrum_hi_radial_s = total_spectrum_hi_radial_s.to(u.K, jybm_to_K)

spec_shape = hi_cube.shape[0]

rot_stack_n = SpectralCube(data=total_spectrum_hi_radial_n.T.reshape((spec_shape, bin_centers.size, 1)),
                           wcs=hi_cube.wcs)
rot_stack_s = SpectralCube(data=total_spectrum_hi_radial_s.T.reshape((spec_shape, bin_centers.size, 1)),
                           wcs=hi_cube.wcs)
rot_stack = SpectralCube(data=total_spectrum_hi_radial.T.reshape((spec_shape, bin_centers.size, 1)),
                         wcs=hi_cube.wcs)

rot_stack_n.write(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_north_{}.fits".format(wstring),
                                              no_check=True), overwrite=True)
rot_stack_s.write(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_south_{}.fits".format(wstring),
                                              no_check=True), overwrite=True)
rot_stack.write(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_{}.fits".format(wstring),
                                            no_check=True), overwrite=True)

np.save(fourteenB_HI_data_wGBT_path("stacked_spectra/radial_stacking_pixelsinbin_{}.npy".format(wstring),
                                    no_check=True),
        num_pixels)
np.save(fourteenB_HI_data_wGBT_path("stacked_spectra/radial_stacking_pixelsinbin_north_{}.npy".format(wstring),
                                    no_check=True),
        num_pixels_n)
np.save(fourteenB_HI_data_wGBT_path("stacked_spectra/radial_stacking_pixelsinbin_south_{}.npy".format(wstring),
                                    no_check=True),
        num_pixels_s)


del hi_cube

hi_cube_cent = SpectralCube.read(fourteenB_wGBT_HI_file_dict["CentSub_Cube"])
hi_mask_cent = fits.open(fourteenB_wGBT_HI_file_dict["CentSub_Mask"])[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

del hi_mask_cent

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

spec_shape = hi_cube_cent.shape[0]
cent_stack_n = SpectralCube(data=total_spectrum_hi_radial_cent_n.T.reshape((spec_shape, bin_centers.size, 1)),
                            wcs=hi_cube_cent.wcs)
cent_stack_s = SpectralCube(data=total_spectrum_hi_radial_cent_s.T.reshape((spec_shape, bin_centers.size, 1)),
                            wcs=hi_cube_cent.wcs)
cent_stack = SpectralCube(data=total_spectrum_hi_radial_cent.T.reshape((spec_shape, bin_centers.size, 1)),
                          wcs=hi_cube_cent.wcs)

cent_stack_n.write(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_north_{}.fits".format(wstring),
                                          no_check=True), overwrite=True)
cent_stack_s.write(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_south_{}.fits".format(wstring),
                                          no_check=True), overwrite=True)
cent_stack.write(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_radial_{}.fits".format(wstring),
                                        no_check=True), overwrite=True)

del hi_cube_cent

hi_cube_peakvel = SpectralCube.read(fourteenB_wGBT_HI_file_dict["PeakSub_Cube"])
hi_mask_peakvel = fits.open(fourteenB_wGBT_HI_file_dict["PeakSub_Mask"])[0]
hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask_peakvel.data > 0)

del hi_mask_peakvel

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

spec_shape = hi_cube_peakvel.shape[0]

peakvel_stack_n = SpectralCube(data=total_spectrum_hi_radial_peakvel_n.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack_s = SpectralCube(data=total_spectrum_hi_radial_peakvel_s.T.reshape((spec_shape, bin_centers.size, 1)),
                               wcs=hi_cube_peakvel.wcs)
peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                             wcs=hi_cube_peakvel.wcs)

peakvel_stack_n.write(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring),
                                                  no_check=True), overwrite=True)
peakvel_stack_s.write(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring),
                                                  no_check=True), overwrite=True)
peakvel_stack.write(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring),
                                                no_check=True), overwrite=True)

del hi_cube_peakvel
