
from spectral_cube import SpectralCube, Projection
import numpy as np
import astropy.units as u
from astropy.io import fits

from spectral_cube.cube_utils import average_beams
from astropy.coordinates import Angle

from paths import (fourteenB_HI_data_path, arecibo_HI_data_path,
                   c_hi_analysispath, paper1_figures_path,
                   data_path, paper1_tables_path)

from constants import (lwidth_name, rotsub_cube_name,
                       rotsub_mask_name, hi_freq,
                       centroidsub_cube_name, centroidsub_mask_name,
                       peakvelsub_cube_name, peakvelssub_mask_name,
                       peaktemps_name)

from plotting_styles import default_figure, onecolumn_figure
from HI_radial_stacking import total_profile


hi_cube = SpectralCube.read(fourteenB_HI_data_path(rotsub_cube_name))
hi_mask = fits.open(fourteenB_HI_data_path(rotsub_mask_name))[0]
hi_cube = hi_cube.with_mask(hi_mask.data > 0)

hi_cube_cent = SpectralCube.read(fourteenB_HI_data_path(centroidsub_cube_name))
hi_mask_cent = fits.open(fourteenB_HI_data_path(centroidsub_mask_name))[0]
hi_cube_cent = hi_cube_cent.with_mask(hi_mask_cent.data > 0)

hi_cube_peakvel = SpectralCube.read(fourteenB_HI_data_path(peakvelsub_cube_name))
hi_mask_peakvel = fits.open(fourteenB_HI_data_path(peakvelssub_mask_name))[0]
hi_cube_peakvel = hi_cube_cent.with_mask(hi_mask_peakvel.data > 0)

hi_peaktemp_hdu = fits.open(fourteenB_HI_data_path(peaktemps_name))[0]
hi_peaktemp = Projection.from_hdu(hi_peaktemp_hdu)

hi_beam = average_beams(hi_cube.beams)

dperc = 5
unit = hi_peaktemp.unit
inneredge = np.nanpercentile(hi_peaktemp, np.arange(0, 101, dperc)[:-1]) * unit
outeredge = np.nanpercentile(hi_peaktemp, np.arange(0, 101, dperc)[1:]) * unit
# Add something small to the 100th percentile so it is used
outeredge[-1] += 1e-3 * unit

total_spectrum_hi_peak = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_hi_peak_cent = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K
total_spectrum_hi_peak_peakvel = \
    np.zeros((inneredge.size, hi_cube.shape[0])) * u.K

for ctr, (p0, p1) in enumerate(zip(inneredge,
                                   outeredge)):

    print("On bin {} to {} K".format(p0, p1))

    mask = np.logical_and(hi_peaktemp >= p0, hi_peaktemp < p1)

    total_spectrum_hi_peak[ctr] = \
        total_profile(hi_cube, mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_peak_cent[ctr] = \
        total_profile(hi_cube_cent, mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))

    total_spectrum_hi_peak_peakvel[ctr] = \
        total_profile(hi_cube_peakvel, mask).to(u.K, equivalencies=hi_beam.jtok_equiv(hi_freq))


# We'll make mock SpectralCubes from these so it's easy to calculate
# moments and such from
rot_stack = SpectralCube(data=total_spectrum_hi_peak.T.reshape((1178, inneredge.size, 1)),
                         wcs=hi_cube.wcs)

cent_stack = SpectralCube(data=total_spectrum_hi_peak_cent.T.reshape((1178, inneredge.size, 1)),
                          wcs=hi_cube.wcs)

peakvel_stack = SpectralCube(data=total_spectrum_hi_peak_peakvel.T.reshape((1178, inneredge.size, 1)),
                             wcs=hi_cube.wcs)

# Now save all of these for future use.
wstring = "{}percentile".format(int(dperc))
rot_stack.write(fourteenB_HI_data_path("stacked_spectra/rotation_stacked_peak_{}.fits".format(wstring),
                                       no_check=True), overwrite=True)

cent_stack.write(fourteenB_HI_data_path("stacked_spectra/centroid_stacked_peak_{}.fits".format(wstring),
                                        no_check=True), overwrite=True)

peakvel_stack.write(fourteenB_HI_data_path("stacked_spectra/peakvel_stacked_peak_{}.fits".format(wstring),
                                           no_check=True), overwrite=True)
