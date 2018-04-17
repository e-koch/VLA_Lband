
'''
Total and radial stacking of the low-res versions of the cubes

Only making peak and centroid stacked outputs.

Also no N/S split.

'''

import astropy.units as u
from spectral_cube import SpectralCube, OneDSpectrum
import numpy as np
from astropy.io import fits
from os.path import join as osjoin
import os
from astropy import log

from cube_analysis.spectral_stacking import radial_stacking

from paths import (iram_co21_14B088_data_path,
                   fourteenB_HI_data_wGBT_path)
from constants import hi_freq
from galaxy_params import gal_feath as gal


smooth_2beam_hi_path = lambda x: osjoin(fourteenB_HI_data_wGBT_path("smooth_2beam"), x)
smooth_2beam_co_path = lambda x: osjoin(iram_co21_14B088_data_path("smooth_2beam"), x)

co_stackpath = lambda x: osjoin(smooth_2beam_co_path("stacked_spectra"), x)
if not os.path.exists(co_stackpath("")):
    os.mkdir(co_stackpath(""))

verbose = False
dr = 500 * u.pc
max_radius = (7.0 * u.kpc).to(u.pc)
wstring = "{0}{1}".format(int(dr.value), dr.unit)
maxrad_string = "{0}{1}".format(int(max_radius.value), max_radius.unit)

# CO stacking first
co_cube = SpectralCube.read(smooth_2beam_co_path("m33.co21_iram.14B-088_HI.38arcsec.centroid_corrected.fits"),
                            memmap=False)
co_cube.allow_huge_operations = True

log.info("Stacking CO with HI centroid. 2 * beam")
bin_centers, total_spectrum_co_radial, num_pixels = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='cube')

spec_shape = co_cube.shape[0]

cent_stack = SpectralCube(data=total_spectrum_co_radial.T.reshape((spec_shape, bin_centers.size, 1)),
                          wcs=co_cube.wcs)
cent_stack.write(co_stackpath("centroid_stacked_radial_{}.fits".format(wstring)),
                 overwrite=True)

# Separately save the number of pixels in each bin
np.save(co_stackpath("radial_stacking_pixelsinbin_{}.npy").format(wstring), num_pixels)

# Save the total profiles over the inner 7 kpc

total_spectrum_co = total_spectrum_co_radial.sum(0)

oned_wcs = co_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co.value, unit=total_spectrum_co.unit,
             wcs=oned_wcs).write(co_stackpath("centroid_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del co_cube

co_cube_peakvel = \
    SpectralCube.read(smooth_2beam_co_path("m33.co21_iram.14B-088_HI.38arcsec.peakvels_corrected.fits"),
                      memmap=False)
co_cube_peakvel.allow_huge_operations = True

log.info("Stacking CO with HI peak vel. 2 * beam")
total_spectrum_co_radial_peakvel = \
    radial_stacking(gal, co_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='cube')[1]

spec_shape = co_cube_peakvel.shape[0]

peakvel_stack = SpectralCube(data=total_spectrum_co_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                             wcs=co_cube_peakvel.wcs)
peakvel_stack.write(co_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)),
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

hi_stackpath = lambda x: osjoin(smooth_2beam_hi_path("stacked_spectra"), x)
if not os.path.exists(hi_stackpath("")):
    os.mkdir(hi_stackpath(""))

hi_cube = SpectralCube.read(smooth_2beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.centroid_corrected.fits"),
                            memmap=True)
# hi_cube.allow_huge_operations = True

hi_mask = fits.open(smooth_2beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec_source_mask.centroid_corrected.fits"))[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

del hi_mask

hi_beam = hi_cube.beam
jybm_to_K = hi_beam.jtok_equiv(hi_freq)

log.info("Stacking HI with HI centroid. 2 * beam")
bin_centers, total_spectrum_hi_radial, num_pixels = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='slice')

total_spectrum_hi_radial = total_spectrum_hi_radial.to(u.K, jybm_to_K)

spec_shape = hi_cube.shape[0]

cent_stack = SpectralCube(data=total_spectrum_hi_radial.T.reshape((spec_shape, bin_centers.size, 1)),
                          wcs=hi_cube.wcs)

cent_stack.write(hi_stackpath("centroid_stacked_radial_{}.fits".format(wstring)),
                 overwrite=True)

np.save(hi_stackpath("radial_stacking_pixelsinbin_{}.npy".format(wstring)),
        num_pixels)

total_spectrum_hi_cent = total_spectrum_hi_radial.sum(0)

# Save each of these
oned_wcs = hi_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_cent.value,
             unit=total_spectrum_hi_cent.unit,
             wcs=oned_wcs).write(hi_stackpath("centroid_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del hi_cube

hi_cube_peakvel = SpectralCube.read(smooth_2beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec.peakvels_corrected.fits"),
    memmap=True)
# hi_cube_peakvel.allow_huge_operations = True

hi_mask_peakvel = fits.open(smooth_2beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.38arcsec_source_mask.peakvels_corrected.fits"))[0]

hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask_peakvel.data > 0)

del hi_mask_peakvel

log.info("Stacking HI with HI peak vel. 2 * beam")
total_spectrum_hi_radial_peakvel = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='slice')[1].to(u.K, jybm_to_K)

spec_shape = hi_cube_peakvel.shape[0]

peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                             wcs=hi_cube_peakvel.wcs)

peakvel_stack.write(hi_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)),
                    overwrite=True)

total_spectrum_hi_peakvel = total_spectrum_hi_radial_peakvel.sum(0)

# Save each of these
oned_wcs = hi_cube_peakvel[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_peakvel.value,
             unit=total_spectrum_hi_peakvel.unit,
             wcs=oned_wcs).write(hi_stackpath("peakvel_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del hi_cube_peakvel


# Onto the 5 * beam versions
log.info("Stacking 5 * beam versions")


smooth_5beam_hi_path = lambda x: osjoin(fourteenB_HI_data_wGBT_path("smooth_5beam"), x)
smooth_5beam_co_path = lambda x: osjoin(iram_co21_14B088_data_path("smooth_5beam"), x)

co_stackpath = lambda x: osjoin(smooth_5beam_co_path("stacked_spectra"), x)
if not os.path.exists(co_stackpath("")):
    os.mkdir(co_stackpath(""))

# CO stacking first
co_cube = SpectralCube.read(smooth_5beam_co_path("m33.co21_iram.14B-088_HI.95arcsec.centroid_corrected.fits"),
                            memmap=False)
co_cube.allow_huge_operations = True

log.info("Stacking CO with HI centroid. 5 * beam")
bin_centers, total_spectrum_co_radial, num_pixels = \
    radial_stacking(gal, co_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='cube')

spec_shape = co_cube.shape[0]

cent_stack = SpectralCube(data=total_spectrum_co_radial.T.reshape((spec_shape, bin_centers.size, 1)),
                          wcs=co_cube.wcs)
cent_stack.write(co_stackpath("centroid_stacked_radial_{}.fits".format(wstring)),
                 overwrite=True)

# Separately save the number of pixels in each bin
np.save(co_stackpath("radial_stacking_pixelsinbin_{}.npy").format(wstring), num_pixels)

# Save the total profiles over the inner 7 kpc

total_spectrum_co = total_spectrum_co_radial.sum(0)

oned_wcs = co_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_co.value, unit=total_spectrum_co.unit,
             wcs=oned_wcs).write(co_stackpath("centroid_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del co_cube

co_cube_peakvel = \
    SpectralCube.read(smooth_5beam_co_path("m33.co21_iram.14B-088_HI.95arcsec.peakvels_corrected.fits"),
                      memmap=False)
co_cube_peakvel.allow_huge_operations = True

log.info("Stacking CO with HI peak vel. 5 * beam")
total_spectrum_co_radial_peakvel = \
    radial_stacking(gal, co_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='cube')[1]

spec_shape = co_cube_peakvel.shape[0]

peakvel_stack = SpectralCube(data=total_spectrum_co_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                             wcs=co_cube_peakvel.wcs)
peakvel_stack.write(co_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)),
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

hi_stackpath = lambda x: osjoin(smooth_5beam_hi_path("stacked_spectra"), x)
if not os.path.exists(hi_stackpath("")):
    os.mkdir(hi_stackpath(""))

hi_cube = SpectralCube.read(smooth_5beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.centroid_corrected.fits"), memmap=True)
# hi_cube.allow_huge_operations = True

hi_mask = fits.open(smooth_5beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.95arcsec_source_mask.centroid_corrected.fits"))[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

del hi_mask

hi_beam = hi_cube.beam
jybm_to_K = hi_beam.jtok_equiv(hi_freq)

log.info("Stacking HI with HI centroid. 5 * beam")
total_spectrum_hi_radial, num_pixels = \
    radial_stacking(gal, hi_cube, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='slice')[1:]

total_spectrum_hi_radial = total_spectrum_hi_radial.to(u.K, jybm_to_K)

spec_shape = hi_cube.shape[0]

cent_stack = SpectralCube(data=total_spectrum_hi_radial.T.reshape((spec_shape, bin_centers.size, 1)),
                          wcs=hi_cube.wcs)

cent_stack.write(hi_stackpath("centroid_stacked_radial_{}.fits".format(wstring)),
                 overwrite=True)

np.save(hi_stackpath("radial_stacking_pixelsinbin_{}.npy".format(wstring)),
        num_pixels)

total_spectrum_hi_cent = total_spectrum_hi_radial.sum(0)

# Save each of these
oned_wcs = hi_cube[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_cent.value,
             unit=total_spectrum_hi_cent.unit,
             wcs=oned_wcs).write(hi_stackpath("centroid_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del hi_cube

hi_cube_peakvel = SpectralCube.read(smooth_5beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.95arcsec.peakvels_corrected.fits"), memmap=True)
# hi_cube_peakvel.allow_huge_operations = True

hi_mask_peakvel = fits.open(smooth_5beam_hi_path("M33_14B-088_HI.clean.image.GBT_feathered.95arcsec_source_mask.peakvels_corrected.fits"))[0]

hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask_peakvel.data > 0)

del hi_mask_peakvel

log.info("Stacking HI with HI peak vel. 5 * beam")
total_spectrum_hi_radial_peakvel = \
    radial_stacking(gal, hi_cube_peakvel, dr=dr,
                    max_radius=max_radius,
                    pa_bounds=None,
                    verbose=verbose,
                    how='slice')[1].to(u.K, jybm_to_K)

spec_shape = hi_cube_peakvel.shape[0]

peakvel_stack = SpectralCube(data=total_spectrum_hi_radial_peakvel.T.reshape((spec_shape, bin_centers.size, 1)),
                             wcs=hi_cube_peakvel.wcs)

peakvel_stack.write(hi_stackpath("peakvel_stacked_radial_{}.fits".format(wstring)),
                    overwrite=True)

total_spectrum_hi_peakvel = total_spectrum_hi_radial_peakvel.sum(0)

# Save each of these
oned_wcs = hi_cube_peakvel[:, 0, 0].wcs
OneDSpectrum(total_spectrum_hi_peakvel.value,
             unit=total_spectrum_hi_peakvel.unit,
             wcs=oned_wcs).write(hi_stackpath("peakvel_stacked_{}.fits".format(maxrad_string)),
                                 overwrite=True)

del hi_cube_peakvel
