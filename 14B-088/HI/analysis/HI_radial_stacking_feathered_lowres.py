
'''
Create stacked profiles of the 38" resolution feathered cube. Only centroid
and peak-velocity.
'''

from spectral_cube import SpectralCube, OneDSpectrum
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import Angle
import os
from os.path import join as osjoin
from astropy import log

from cube_analysis.spectral_stacking import total_profile, radial_stacking

from paths import fourteenB_HI_data_wGBT_path

from constants import hi_freq

from galaxy_params import gal

smooth_2beam_folder = lambda x: osjoin(fourteenB_HI_data_wGBT_path("smooth_2beam"), x)

hi_cube = SpectralCube.read(smooth_2beam_folder("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.centroid_corrected.fits"),
                            memmap=True)
hi_beam = hi_cube.beam
jybm_to_K = hi_beam.jtok_equiv(hi_freq)

hi_radius = gal.radius(header=hi_cube.header)
# hi_pas = gal.position_angles(header=hi_cube.header)

verbose = False

dr = 200 * u.pc
max_radius = (10.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)
maxrad_string = "{0}{1}".format(int(max_radius.value), max_radius.unit)

nbins = np.int(np.floor(max_radius / dr))
inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)
# Make the N and S masks
# pa_bounds_n = Angle([0.5 * np.pi * u.rad, -0.5 * np.pi * u.rad])

# pa_bounds_s = Angle([-0.5 * np.pi * u.rad, 0.5 * np.pi * u.rad])

hi_mask = fits.open(smooth_2beam_folder("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec_source_mask.centroid_corrected.fits"))[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

del hi_mask

cent_shape = hi_cube.shape[0]

log.info("Radial stacking centroid")
bin_centers, total_spectrum_hi_radial, num_pixels = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='slice')

total_spectrum_hi_radial = total_spectrum_hi_radial.to(u.K, jybm_to_K)

np.save(smooth_2beam_folder("stacked_spectra/radial_stacking_pixelsinbin_{}.npy").format(wstring), num_pixels)

# total_spectrum_hi_radial_n = \
#     radial_stacking(gal, hi_cube, dr=dr,
#                     max_radius=max_radius,
#                     pa_bounds=pa_bounds_n,
#                     verbose=verbose,
#                     how='slice')[1].to(u.K, jybm_to_K)

# total_spectrum_hi_radial_s = \
#     radial_stacking(gal, hi_cube, dr=dr,
#                     max_radius=max_radius,
#                     pa_bounds=pa_bounds_s,
#                     verbose=verbose,
#                     how='slice')[1].to(u.K, jybm_to_K)

# Now save all of these for future use.
stacked_folder = smooth_2beam_folder("stacked_spectra")
if not os.path.exists(stacked_folder):
    os.mkdir(stacked_folder)

# Radial stacks
# cent_stack_n = SpectralCube(data=total_spectrum_hi_radial_n.T.reshape((cent_shape, inneredge.size, 1)),
#                            wcs=hi_cube.wcs)
# cent_stack_s = SpectralCube(data=total_spectrum_hi_radial_s.T.reshape((cent_shape, inneredge.size, 1)),
#                            wcs=hi_cube.wcs)
cent_stack = SpectralCube(data=total_spectrum_hi_radial.T.reshape((cent_shape, inneredge.size, 1)),
                          wcs=hi_cube.wcs)

# rot_stack_n.write(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_north_{}.fits".format(wstring),
#                                          no_check=True), overwrite=True)
# rot_stack_s.write(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_radial_south_{}.fits".format(wstring),
#                                          no_check=True), overwrite=True)
cent_stack.write(smooth_2beam_folder("stacked_spectra/centroid_stacked_radial_{}.fits".format(wstring)),
                 overwrite=True)

# Total over all spectra
total_spectrum_hi = total_spectrum_hi_radial.sum(0)

oned_wcs = hi_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi.value, unit=total_spectrum_hi.unit,
             wcs=oned_wcs).write(smooth_2beam_folder("stacked_spectra/centroid_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

# Now delete that cube to save memory before loading the next
del hi_cube

hi_cube_peakvel = SpectralCube.read(smooth_2beam_folder("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.peakvels_corrected.fits"),
                                    memmap=True)
hi_mask_peakvel = fits.open(smooth_2beam_folder("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec_source_mask.peakvels_corrected.fits"))[0]

hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask_peakvel.data > 0)

del hi_mask_peakvel

peak_shape = hi_cube_peakvel.shape[0]

log.info("Radial stacking peak vel.")

total_spectrum_hi_radial_peakvel = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='slice')[1].to(u.K, jybm_to_K)

# total_spectrum_hi_radial_peakvel_n = \
#     radial_stacking(gal, hi_cube_peakvel, dr=dr,
#                     max_radius=max_radius,
#                     pa_bounds=pa_bounds_n,
#                     verbose=verbose,
#                     how='slice')[1].to(u.K, jybm_to_K)

# total_spectrum_hi_radial_peakvel_s = \
#     radial_stacking(gal, hi_cube_peakvel, dr=dr,
#                     max_radius=max_radius,
#                     pa_bounds=pa_bounds_s,
#                     verbose=verbose,
#                     how='slice')[1].to(u.K, jybm_to_K)

# peakvel_stack_n = SpectralCube(data=total_spectrum_hi_radial_peakvel_n.T.reshape((peak_shape, inneredge.size, 1)),
#                                wcs=hi_cube_peakvel.wcs)
# peakvel_stack_s = SpectralCube(data=total_spectrum_hi_radial_peakvel_s.T.reshape((peak_shape, inneredge.size, 1)),
#                                wcs=hi_cube_peakvel.wcs)
peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((peak_shape, inneredge.size, 1)),
                             wcs=hi_cube_peakvel.wcs)

# peakvel_stack_n.write(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_north_{}.fits".format(wstring),
#                                                   no_check=True), overwrite=True)
# peakvel_stack_s.write(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_radial_south_{}.fits".format(wstring),
#                                                   no_check=True), overwrite=True)
peakvel_stack.write(smooth_2beam_folder("stacked_spectra/peakvel_stacked_radial_{}.fits".format(wstring)),
                    overwrite=True)

total_spectrum_hi_peakvel = total_spectrum_hi_radial_peakvel.sum(0)

oned_wcs = hi_cube_peakvel[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_peakvel.value, unit=total_spectrum_hi_peakvel.unit,
             wcs=oned_wcs).write(smooth_2beam_folder("stacked_spectra/peakvel_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del hi_cube_peakvel
