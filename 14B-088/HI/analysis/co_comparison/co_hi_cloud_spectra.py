
'''
Use the clean cloud sample from cloud_catalog.py to find the CO and HI spectra
over each cloud.
'''

from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.console import ProgressBar
import astropy.units as u
from scipy import ndimage as nd
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.cube_utils import average_beams
from reproject import reproject_interp
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import seaborn as sb
from signal_id import Noise

from cube_analysis.spectra_shifter import cube_shifter

from analysis.paths import (fourteenB_HI_data_path, iram_co21_data_path,
                            paper1_figures_path)
from analysis.constants import hi_freq, cube_name, moment1_name
from analysis.galaxy_params import gal
from analysis.plotting_styles import onecolumn_figure, align_yaxis


hi_cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))
hi_beam = average_beams(hi_cube.beams)
hi_mom1 = fits.open(fourteenB_HI_data_path(moment1_name))[0]


# cloud_mask = \
#     fits.open(iram_co21_data_path("m33.co21_new_assign_cprops_cleansample.fits"))[0]
# Improved cloud mask from Braine & Corbelli
cloud_mask_hdu = fits.open(iram_co21_data_path("asgn.fits"))[0]
# The cloud edges can be quite strict. Spectrally expand the mask one pixel
# in each direction
cloud_mask = nd.binary_dilation(cloud_mask_hdu.data > 0,
                                np.array([[[1, 1, 1]]]).T)

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
cube = cube.with_mask(cloud_mask)

# Make a SNR cube
# noise = Noise(cube)
# noise.estimate_noise(spectral_flat=True)
# noise.get_scale_cube()

# snr = SpectralCube(data=noise.snr.copy(), wcs=cube.wcs)


# Create a CO centroid map where the GMC mask is valid.
def peak_velocity(y, x, cube):
    argmax = np.nanargmax(cube[:, y, x].value)
    return cube.spectral_axis[argmax]


peak_intens = cube.max(0)
valid_mask = cloud_mask.sum(0) > 0

cube_specinterp = cube.spectral_interpolate(hi_cube.spectral_axis)
mom1 = cube_specinterp.moment1()
# mom1[~peak_mask] = np.NaN
peak_vels_arr = np.zeros_like(mom1) * np.NaN
pbar = ProgressBar(len(zip(*np.where(valid_mask))))
for y, x in zip(*np.where(valid_mask)):
    peak_vels_arr[y, x] = peak_velocity(y, x, cube_specinterp)
    pbar.update()

# Moment 1 and peak vels need to be somewhat close to each other
# max_diff = 20 * u.km / u.s
# bad_peaks = np.abs(peak_vels_arr - mom1) > max_diff
# peak_vels_arr[bad_peaks] = np.NaN

peak_vels = Projection(peak_vels_arr, unit=mom1.unit, wcs=mom1.wcs)
peak_vels.write(iram_co21_data_path("m33.co21_iram.masked.peakvel.fits", no_check=True))

mom1_reproj_arr = reproject_interp(mom1.hdu,
                                   hi_cube.wcs.celestial,
                                   shape_out=hi_cube.shape[1:])[0]
mom1_reproj = Projection(mom1_reproj_arr, unit=mom1.unit,
                         wcs=hi_cube.wcs.celestial)
peak_vels_reproj_arr = \
    reproject_interp(peak_vels.hdu,
                     hi_cube.wcs.celestial,
                     shape_out=hi_cube.shape[1:])[0]
peak_vels_reproj = Projection(peak_vels_reproj_arr, unit=mom1.unit,
                              wcs=hi_cube.wcs.celestial)
peak_vels_reproj.write(iram_co21_data_path("m33.co21_iram.masked.peakvel.hireprojection.fits", no_check=True))

# hi_mom1_reproj_arr = reproject_interp(hi_mom1, cube.wcs.celestial,
#                                       shape_out=mom1.shape)[0]
# hi_mom1_reproj = Projection(hi_mom1_reproj_arr, unit=mom1.unit,
#                             wcs=cube.wcs.celestial)

# vels_co = hi_mom1_reproj
# vels_hi = Projection(hi_mom1.data, header=hi_mom1.header, unit=u.m / u.s)

vels_co = peak_vels
vels_hi = peak_vels_reproj

# vels_co = mom1
# vels_hi = mom1_reproj

# Loop through the clouds, making the total/average profiles.
num_clouds = int(cloud_mask_hdu.data.max())

cloud_total_specs_co = np.zeros((num_clouds, cube.shape[0]))
cloud_total_specs_hi = np.zeros((num_clouds, hi_cube.shape[0]))

cloud_avg_specs_co = np.zeros((num_clouds, cube.shape[0]))
cloud_avg_specs_hi = np.zeros((num_clouds, hi_cube.shape[0]))

# Plot for each cloud
verbose = False

# pool = Pool(6)
pool = None

for i in ProgressBar(num_clouds):

    # Find spatial extent of the cloud
    plane = (cloud_mask_hdu.data == i + 1).sum(0) > 0
    xy_posns = np.where(plane)

    co_spec = np.zeros((cube.shape[0],))
    co_count = np.zeros((cube.shape[0],))

    # Shift wrt to CO peak or centroid.
    co_shifted = cube_shifter(cube, vels_co, gal.vsys, xy_posns=xy_posns,
                              pool=pool, return_spectra=True,)

    for spec in co_shifted[0]:
        finites = np.isfinite(spec)
        co_spec[finites] += spec.value[finites]
        co_count += finites

    cloud_total_specs_co[i] = co_spec
    cloud_avg_specs_co[i] = co_spec / co_count

    # Now regrid plane onto the HI and do the same.
    plane_reproj = reproject_interp((plane, WCS(cloud_mask_hdu.header).celestial),
                                    hi_cube.wcs.celestial,
                                    shape_out=hi_cube.shape[1:])[0] > 0
    xy_posns_hi = np.where(plane_reproj)

    hi_spec = np.zeros((hi_cube.shape[0],))
    hi_count = np.zeros((hi_cube.shape[0],))

    hi_shifted = cube_shifter(hi_cube, vels_hi, gal.vsys,
                              xy_posns=xy_posns_hi,
                              pool=pool, return_spectra=True,)

    for spec in hi_shifted[0]:
        # if np.nanmax(spec.value) < 0.0066:
        #     continue
        finites = np.isfinite(spec) & (spec > 0)
        hi_spec[finites] += spec.value[finites]
        hi_count += finites

    cloud_total_specs_hi[i] = (hi_spec * hi_beam.jtok(hi_freq).value)
    cloud_avg_specs_hi[i] = (hi_spec * hi_beam.jtok(hi_freq).value) / hi_count

    if verbose:

        onecolumn_figure()
        sb.set_context("paper", rc={"font.family": "serif",
                                    "figure.figsize": np.array([7.325, 5.9375])},
                       font_scale=1.5)
        sb.set_palette("afmhot")

        fig, ax = plt.subplots(3, 1,
                               sharex=True,
                               sharey=False, num=1)

        plt.subplots_adjust(hspace=0.1,
                            wspace=0.1)

        fig.text(0.5, 0.02, 'Velocity (km/s)', ha='center')

        ax[0].plot(hi_cube.spectral_axis.to(u.km / u.s).value,
                   cloud_avg_specs_hi[i], 'b',
                   drawstyle='steps-mid', label='HI')
        # For the legend
        ax[0].plot([], [], 'r--', label="CO(2-1)")
        ax2 = ax[0].twinx()
        ax2.plot(cube.spectral_axis.to(u.km / u.s).value,
                 cloud_avg_specs_co[i] * 1000.,
                 'r--',
                 drawstyle='steps-mid', label='CO(2-1)')
        align_yaxis(ax[0], 0, ax2, 0)
        # ax[0].set_ylim([0, 1])
        ax[0].set_xlim([gal.vsys.to(u.km / u.s).value - 40,
                        gal.vsys.to(u.km / u.s).value + 40])
        ax[0].set_ylabel("HI Intensity (K)")
        ax2.set_ylabel("CO Intensity (mK)")
        ax[0].grid()

        ax[0].legend(frameon=True, loc='upper right')

        for spec in hi_shifted[0]:
            if not np.isnan(spec).all():
                ax[1].plot(hi_cube.spectral_axis.to(u.km / u.s).value,
                           spec.value * hi_beam.jtok(hi_freq).value, '--',
                           alpha=0.5, drawstyle='steps-mid')
        ax[1].set_ylabel("HI Intensity (K)")
        ax[1].grid()

        for spec in co_shifted[0]:
            if not np.isnan(spec).all():
                ax[2].plot(cube.spectral_axis.to(u.km / u.s).value,
                           spec.to(u.mK).value, '--',
                           alpha=0.5, drawstyle='steps-mid')
        ax[2].set_ylabel("CO Intensity (mK)")
        ax[2].grid()

        # plt.draw()
        # raw_input("?")
        # plt.clf()
        fig.savefig(paper1_figures_path("GMC_CO_peakshift/GMC_CO_peakshift"
                                        "_{}.png".format(i + 1)))
        fig.savefig(paper1_figures_path("GMC_CO_peakshift/GMC_CO_peakshift"
                                        "_{}.pdf".format(i + 1)))

        plt.clf()
plt.close()

# Save the stacked spectra
np.savetxt("GMC_stackedspectra_peakvels_CO.txt", cloud_avg_specs_co)
np.savetxt("GMC_stackedspectra_peakvels_HI.txt", cloud_avg_specs_hi)
