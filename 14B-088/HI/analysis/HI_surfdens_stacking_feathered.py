
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

'''
Stack the peak vel subtracted cube in intervals of surface density
'''

hi_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict["PeakSub_Cube"])
hi_mask = fits.open(fourteenB_wGBT_HI_file_dict["PeakSub_Mask"])[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

del hi_mask

hi_mom0_hdu = fits.open(fourteenB_HI_data_wGBT_path("M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.moment0_Kkms.fits"))[0]
hi_mom0 = Projection.from_hdu(hi_mom0_hdu)

hi_beam = hi_cube.beam

dperc = 5
unit = hi_mom0.unit
inneredge = np.nanpercentile(hi_mom0, np.arange(0, 101, dperc)[:-1]) * unit
outeredge = np.nanpercentile(hi_mom0, np.arange(0, 101, dperc)[1:]) * unit
# Add something small to the 100th percentile so it is used
outeredge[-1] += 1e-3 * unit

total_spectrum_hi_peak = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K

num_pix = np.zeros_like(inneredge.value)

log.info("Running peak sub")
for ctr, (p0, p1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {} K km/s".format(p0, p1))

    mask = np.logical_and(hi_mom0 >= p0, hi_mom0 < p1)

    num_pix[ctr] = float(mask.sum())

    total_spectrum_hi_peak[ctr] = \
        total_profile(hi_cube, mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

mask = BooleanArrayMask(np.ones((hi_cube.shape[0], inneredge.size, 1), dtype=bool),
                        hi_cube.wcs)
peakvel_stack = SpectralCube(data=total_spectrum_hi_peak.T.reshape((hi_cube.shape[0], inneredge.size, 1)),
                             wcs=hi_cube.wcs, mask=mask)

stacked_folder = fourteenB_HI_data_wGBT_path("stacked_spectra", no_check=True)
if not os.path.exists(stacked_folder):
    os.mkdir(stacked_folder)

wstring = "{}percentile".format(int(dperc))
peakvel_stack.write(fourteenB_HI_data_wGBT_path("stacked_spectra/peakvel_stacked_surfdens_{}.fits".format(wstring),
                                            no_check=True), overwrite=True)

# Save the number of pixels in each bin
np.save(fourteenB_HI_data_wGBT_path("stacked_spectra/surfdens_stacking_{}_num_pix.npy".format(wstring),
                                    no_check=True), num_pix)

del hi_cube
