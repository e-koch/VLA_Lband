
from spectral_cube import SpectralCube, Projection, BooleanArrayMask
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy import log

import os

from cube_analysis.spectral_stacking import total_profile

from paths import (fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict,
                   allfigs_path, alltables_path)

from constants import hi_freq


hi_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["RotSub_Cube"])
hi_mask = fits.open(fourteenB_wGBT_HI_file_dict["RotSub_Mask"])[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

del hi_mask

hi_peaktemp_hdu = fits.open(fourteenB_wGBT_HI_file_dict["PeakTemp"])[0]
hi_peaktemp = Projection.from_hdu(hi_peaktemp_hdu)

hi_beam = hi_cube.beam

dperc = 5
unit = hi_peaktemp.unit
inneredge = np.nanpercentile(hi_peaktemp, np.arange(0, 101, dperc)[:-1]) * unit
outeredge = np.nanpercentile(hi_peaktemp, np.arange(0, 101, dperc)[1:]) * unit
# Add something small to the 100th percentile so it is used
outeredge[-1] += 1e-3 * unit

total_spectrum_hi_peak = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K

num_pix = np.zeros_like(inneredge.value)

log.info("Running rot sub")
for ctr, (p0, p1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {} K".format(p0, p1))

    mask = np.logical_and(hi_peaktemp >= p0, hi_peaktemp < p1)

    num_pix[ctr] = float(mask.sum())

    total_spectrum_hi_peak[ctr] = \
        total_profile(hi_cube, mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

mask = BooleanArrayMask(np.ones((hi_cube.shape[0], inneredge.size, 1), dtype=bool),
                        hi_cube.wcs)
rot_stack = SpectralCube(data=total_spectrum_hi_peak.T.reshape((hi_cube.shape[0], inneredge.size, 1)),
                         wcs=hi_cube.wcs, mask=mask)

stacked_folder = fourteenB_HI_data_wGBT_path("stacked_spectra", no_check=True)
if not os.path.exists(stacked_folder):
    os.mkdir(stacked_folder)

wstring = "{}percentile".format(int(dperc))
rot_stack.write(fourteenB_HI_data_wGBT_path("stacked_spectra/rotation_stacked_peak_{}.fits".format(wstring),
                                            no_check=True), overwrite=True)

# Save the number of pixels in each bin
np.save(fourteenB_HI_data_wGBT_path("stacked_spectra/peak_stacking_{}_num_pix.npy".format(wstring),
                                    no_check=True), num_pix)

del hi_cube

hi_cube_cent = SpectralCube.read(fourteenB_wGBT_HI_file_dict["CentSub_Cube"])
hi_mask_cent = fits.open(fourteenB_wGBT_HI_file_dict["CentSub_Mask"])[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

total_spectrum_hi_peak_cent = \
    np.zeros((inneredge.size, hi_cube_cent.shape[0])) * u.K

del hi_mask_cent

log.info("Running cent sub")
for ctr, (p0, p1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {} K".format(p0, p1))

    mask = np.logical_and(hi_peaktemp >= p0, hi_peaktemp < p1)

    total_spectrum_hi_peak_cent[ctr] = \
        total_profile(hi_cube_cent, mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

mask = BooleanArrayMask(np.ones((hi_cube_cent.shape[0], inneredge.size, 1), dtype=bool),
                        hi_cube_cent.wcs)
cent_stack = SpectralCube(data=total_spectrum_hi_peak_cent.T.reshape((hi_cube_cent.shape[0], inneredge.size, 1)),
                          wcs=hi_cube_cent.wcs, mask=mask)

cent_stack.write(fourteenB_HI_data_wGBT_path("stacked_spectra/centroid_stacked_peak_{}.fits".format(wstring),
                                             no_check=True), overwrite=True)

del hi_cube_cent

hi_cube_peakvel = SpectralCube.read(fourteenB_wGBT_HI_file_dict["PeakSub_Cube"])
hi_mask_peakvel = fits.open(fourteenB_wGBT_HI_file_dict["PeakSub_Mask"])[0]
hi_cube_peakvel = hi_cube_peakvel.with_mask(hi_mask_peakvel.data > 0)

total_spectrum_hi_peak_peakvel = \
    np.zeros((inneredge.size, hi_cube_peakvel.shape[0])) * u.K

del hi_mask_peakvel

log.info("Running peak sub")
for ctr, (p0, p1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {} K".format(p0, p1))

    mask = np.logical_and(hi_peaktemp >= p0, hi_peaktemp < p1)

    total_spectrum_hi_peak_peakvel[ctr] = \
        total_profile(hi_cube_peakvel, mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

mask = BooleanArrayMask(np.ones((hi_cube_peakvel.shape[0], inneredge.size, 1), dtype=bool),
                        hi_cube_peakvel.wcs)
peakvel_stack = SpectralCube(data=total_spectrum_hi_peak_peakvel.T.reshape((hi_cube_peakvel.shape[0], inneredge.size, 1)),
                             wcs=hi_cube_peakvel.wcs, mask=mask)

peakvel_stack.write(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_peak_{}.fits".format(wstring),
                                                no_check=True), overwrite=True)

del hi_cube_peakvel
