
'''
Use the clean cloud sample from cloud_catalog.py to find the CO and HI spectra
over each cloud.
'''

from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.console import ProgressBar
import astropy.units as u
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.cube_utils import average_beams
from reproject import reproject_interp
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool


from analysis.paths import (fourteenB_HI_data_path, iram_co21_data_path,
                            paper1_figures_path)
from analysis.constants import rotsub_cube_name, hi_freq, cube_name, moment1_name
from analysis.galaxy_params import gal
from analysis.spectra_shifter import cube_shifter


# Set the same origin for double axis plots
def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1 - y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny + dy, maxy + dy)


hi_cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))
hi_beam = average_beams(hi_cube.beams)
hi_mom1 = fits.open(fourteenB_HI_data_path(moment1_name))[0]

tab = Table.read(iram_co21_data_path("m33.co21_new_props_clfind_cleansample.fits"))

cloud_mask = \
    fits.open(iram_co21_data_path("m33.co21_new_assign_cprops_cleansample.fits"))[0]

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
del cube._header[""]


# Create a CO centroid map where the GMC mask is valid.
def peak_velocity(y, x, cube):
    argmax = np.argmax(cube[:, y, x].value)
    return cube.spectral_axis[argmax]


peak_intens = cube.max(0)
peak_mask = peak_intens > 5 * 22. * u.mK
valid_mask = cloud_mask.data.sum(0) > 0

cube_specinterp = cube.spectral_interpolate(hi_cube.spectral_axis)
mom1 = cube_specinterp.with_mask(valid_mask).moment1()
mom1[~peak_mask] = np.NaN
# mom1[~peak_mask] = np.NaN
peak_vels_arr = np.zeros_like(mom1)
for y, x in zip(*np.where(valid_mask)):
    peak_vels_arr[y, x] = peak_velocity(y, x, cube_specinterp)
peak_vels_arr[peak_vels_arr == 0.0] = np.NaN
peak_vels_arr[~peak_mask] = np.NaN
# Moment 1 and peak vels need to be somewhat close to each other
# max_diff = 20 * u.km / u.s
# bad_peaks = np.abs(peak_vels_arr - mom1) > max_diff
# peak_vels_arr[bad_peaks] = np.NaN

peak_vels = Projection(peak_vels_arr, unit=mom1.unit, wcs=mom1.wcs)

# mom1_reproj_arr = reproject_interp(mom1.hdu,
#                                    hi_cube.wcs.celestial,
#                                    shape_out=hi_cube.shape[1:])[0]
# mom1_reproj = Projection(mom1_reproj_arr, unit=mom1.unit,
#                          wcs=hi_cube.wcs.celestial)
peak_vels_reproj_arr = \
    reproject_interp(peak_vels.hdu,
                     hi_cube.wcs.celestial,
                     shape_out=hi_cube.shape[1:])[0]
peak_vels_reproj = Projection(peak_vels_reproj_arr, unit=mom1.unit,
                              wcs=hi_cube.wcs.celestial)

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

# Reproject onto HI grid

# Loop through the clouds, making the total/average profiles.
num_clouds = int(cloud_mask.data.max())

cloud_total_specs_co = np.zeros((num_clouds, cube.shape[0]))
cloud_total_specs_hi = np.zeros((num_clouds, hi_cube.shape[0]))

cloud_avg_specs_co = np.zeros((num_clouds, cube.shape[0]))
cloud_avg_specs_hi = np.zeros((num_clouds, hi_cube.shape[0]))

# Plot for each cloud
verbose = True

# pool = Pool(6)
pool = None

for i in ProgressBar(num_clouds):

    # Find spatial extent of the cloud
    plane = (cloud_mask.data == i + 1).sum(0) > 0
    xy_posns = np.where(plane)

    co_spec = np.zeros((cube.shape[0],))
    co_count = np.zeros((cube.shape[0],))

    # Shift wrt to CO centroids.
    co_shifted = cube_shifter(cube, vels_co, gal.vsys, xy_posns=xy_posns,
                              pool=pool, return_spectra=True,)

    for spec in co_shifted[0]:
        finites = np.isfinite(spec)
        co_spec[finites] += spec.value[finites]
        co_count += finites

    cloud_total_specs_co[i] = co_spec
    cloud_avg_specs_co[i] = co_spec / co_count

    # Now regrid plane onto the HI and do the same.
    plane_reproj = reproject_interp((plane, WCS(cloud_mask.header).celestial),
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
        fig, ax = plt.subplots(3, 1,
                               sharex=True,
                               sharey=False, num=1)

        plt.subplots_adjust(hspace=0.1,
                            wspace=0.1)

        fig.text(0.5, 0.02, 'Velocity (km/s)', ha='center')

        ax[0].plot(hi_cube.spectral_axis.value,
                   cloud_avg_specs_hi[i], 'b',
                   drawstyle='steps-mid', label='HI')
        # For the legend
        ax[0].plot([], [], 'r--', label="CO(2-1)")
        ax2 = ax[0].twinx()
        ax2.plot(cube.spectral_axis.value,
                 cloud_avg_specs_co[i],
                 'r--',
                 drawstyle='steps-mid', label='CO(2-1)')
        align_yaxis(ax[0], 0, ax2, 0)
        # ax[0].set_ylim([0, 1])
        ax[0].set_xlim([gal.vsys.to(u.m / u.s).value - 40000,
                        gal.vsys.to(u.m / u.s).value + 40000])
        ax[0].set_ylabel("HI Intensity (K)")
        ax2.set_ylabel("CO Intensity (K)")
        ax[0].grid()

        ax[0].legend()

        for spec in hi_shifted[0]:
            if not np.isnan(spec).all():
                ax[1].plot(hi_cube.spectral_axis.value,
                           spec.value * hi_beam.jtok(hi_freq).value, '--',
                           alpha=0.5, drawstyle='steps-mid')
        ax[1].set_ylabel("HI Intensity (K)")
        ax[1].grid()

        for spec in co_shifted[0]:
            if not np.isnan(spec).all():
                ax[2].plot(cube.spectral_axis.value, spec.value, '--',
                           alpha=0.5, drawstyle='steps-mid')
        ax[2].set_ylabel("CO Intensity (K)")
        ax[2].grid()

        # plt.draw()

        # raw_input("?")
        fig.savefig(paper1_figures_path("GMC_CO_peakshift/GMC_CO_peakshift"
                                        "_{}.png".format(i + 1)))

        plt.clf()
plt.close()
