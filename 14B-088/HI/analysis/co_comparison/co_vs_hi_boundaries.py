from basics import BubbleFinder2D
from basics.utils import sig_clip
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.cube_utils import average_beams
import astropy.units as u
import matplotlib.pyplot as p
import scipy.ndimage as nd
from scipy.stats import binned_statistic
import numpy as np
from skimage.morphology import medial_axis
from reproject import reproject_interp
from skimage.segmentation import find_boundaries
from astropy.utils.console import ProgressBar

from analysis.paths import (fourteenB_HI_data_path, iram_co21_data_path,
                            paper1_figures_path)
from analysis.constants import rotsub_cube_name, hi_freq
from analysis.galaxy_params import gal

'''
Calculate the intensities of HI and CO as a function of distance from the
edges of the adap. thresh mask.

To avoid a whole cube regrid, I use the rotation subtracted cubes and extract
zeroth moments from spectral slabs.

'''

np.random.seed(34678953)

# Plot a bunch
verbose = False
# slicer = (slice(825, 1033), slice(360, 692))
slicer = (slice(None), slice(None))

# Load in the rotation subtracted cubes
hi_cube = SpectralCube.read(fourteenB_HI_data_path(rotsub_cube_name))
co_cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.rotsub.fits"))

start_vel = - 30 * u.km / u.s
end_vel = 30 * u.km / u.s
del_vel = 3. * u.km / u.s

vels = np.arange(start_vel.value, end_vel.value + del_vel.value,
                 del_vel.value) * u.km / u.s

# Get the radius array so we can cut to where the CO data is valid
radii = gal.radius(header=hi_cube[0].header)
max_radius = 6.0 * u.kpc

all_dists = []
all_radii = []
all_vals_hi = []
all_vals_co = []
edge_masks = []
skeletons = []

hi_beam = average_beams(hi_cube.beams)

# Estimate the noise level in an equivalent slab
hi_mom0 = hi_cube.spectral_slab(-180 * u.km / u.s, -183 * u.km / u.s).moment0()
sigma = sig_clip(hi_mom0.value, nsig=10) * \
    hi_beam.jtok(hi_freq).value

for i, (end, start) in enumerate(ProgressBar(zip(vels[1:], vels[:-1]))):

    hi_slab = hi_cube.spectral_slab(start, end)
    hi_mom0 = hi_slab.moment0() * hi_beam.jtok(hi_freq) / u.Jy

    # Make the CO slab, then reproject onto the HI grid
    co_mom0 = co_cube.spectral_slab(start, end).moment0()
    co_mom0_reproj = reproject_interp(co_mom0.hdu, hi_mom0.header)[0]
    co_mom0 = Projection(co_mom0_reproj, wcs=hi_mom0.wcs)

    # Need a mask from the HI
    # Adjust the sigma in a single channel to the moment0 in the slab
    # sigma = 0.00152659 * hi_slab.shape[0] * \
    #     np.abs((hi_slab.spectral_axis[1] - hi_slab.spectral_axis[0]).value)

    bub = BubbleFinder2D(hi_mom0, auto_cut=False, sigma=sigma)
    bub.create_mask(bkg_nsig=30, region_min_nsig=60, mask_clear_border=False)

    # skeleton = medial_axis(~bub.mask)
    # skeletons.append(skeleton)

    edge_mask = find_boundaries(bub.mask, connectivity=2, mode='outer')
    hole_mask = bub.mask.copy()
    # Now apply a radial boundary to the edge mask where the CO data is valid
    # This is the same cut-off used to define the valid clouds
    radial_cut = radii <= max_radius
    edge_mask *= radial_cut
    edge_masks.append(edge_mask)

    dist_trans = nd.distance_transform_edt(~edge_mask)
    # Assign negative values to regions within holes.
    dist_trans[hole_mask] = -dist_trans[hole_mask]

    # hist = p.hist(co_mom0.value[np.isfinite(co_mom0.value)], bins=100)
    # p.draw()
    # raw_input("?")
    # p.clf()

    # print(np.nanmin(co_mom0.value))

    # print(np.nanmax(co_mom0.value[radial_cut]))
    # print(np.nanmin(co_mom0.value[radial_cut]))

    all_dists.extend(list(dist_trans[radial_cut]))
    all_radii.extend(list(radii.value[radial_cut]))
    all_vals_hi.extend(list(hi_mom0.value[radial_cut] / 1000.))
    all_vals_co.extend(list(co_mom0.value[radial_cut] / 1000.))

    if verbose:
        print("Velocities: {} to {}".format(start, end))
        ax = p.subplot(121, projection=hi_mom0.wcs)
        p.imshow(np.arctan(hi_mom0[slicer].value /
                           np.nanpercentile(hi_mom0.value, 85)),
                 origin='lower')
        # p.contour(skeleton, colors='b')
        p.contour(edge_mask[slicer], colors='g')

        ax2 = p.subplot(122, projection=hi_mom0.wcs)
        # p.imshow(hole_mask, origin='lower')
        p.imshow(np.arctan(co_mom0.value[slicer] /
                           np.nanpercentile(co_mom0.value, 95)),
                 origin='lower')
        # p.imshow(co_mom0.value, origin='lower')
        # p.contour(skeleton, colors='b')
        p.contour(edge_mask[slicer], colors='g')
        p.draw()
        lat = ax2.coords[1]
        lat.set_ticklabel_visible(False)

        raw_input("Next plot?")
        p.clf()

all_dists = np.array(all_dists)
all_radii = np.array(all_radii)
all_vals_hi = np.array(all_vals_hi)
all_vals_co = np.array(all_vals_co)

# Now bin all of the distances against the HI and CO intensities.
bins = np.arange(-30, 30, 1)
hi_vals, bin_edges, bin_num = \
    binned_statistic(all_dists, all_vals_hi,
                     bins=bins,
                     statistic=np.nanmean)

co_vals, bin_edges, bin_num = \
    binned_statistic(all_dists, all_vals_co,
                     bins=bins,
                     statistic=np.nanmean)

binned_elements = \
    binned_statistic(all_dists, np.ones_like(all_dists), bins=bins,
                     statistic=np.sum)[0]

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width / 2

# Let's bootstrap to get errors in the distance bins
niters = 100
hi_samps = np.zeros((niters, len(bin_centers)))
co_samps = np.zeros((niters, len(bin_centers)))
print("Bootstrapping")
for i in ProgressBar(niters):
    hi_samps[i] = \
        binned_statistic(all_dists, np.random.permutation(all_vals_hi),
                         bins=bins,
                         statistic=np.nanmean)[0]

    co_samps[i] = \
        binned_statistic(all_dists, np.random.permutation(all_vals_co),
                         bins=bins,
                         statistic=np.nanmean)[0]

# Take the stds in the distribution for each bin
hi_errs = np.nanstd(hi_samps, axis=0)
co_errs = np.nanstd(co_samps, axis=0)

# Convert the bin_centers to pc
pixscale = \
    hi_mom0.header['CDELT2'] * (np.pi / 180.) * gal.distance.to(u.pc).value

bin_centers *= pixscale

p.errorbar(bin_centers, hi_vals / np.nanmax(hi_vals),
           yerr=hi_errs / np.nanmax(hi_vals), fmt="D-",
           color="b", label="HI")
p.errorbar(bin_centers, co_vals / np.nanmax(co_vals),
           yerr=co_errs / np.nanmax(co_vals), fmt="o-",
           color="r", label="CO(2-1)")

# p.plot(bin_centers, hi_vals / np.nanmax(hi_vals), 'bD-',
#        label="HI")
# p.plot(bin_centers, co_vals / np.nanmax(co_vals), 'ro-',
#        label="CO(2-1)")
# p.xlim([0.0, 200])
p.ylim([0.0, 1.1])
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Normalized Intensity")
p.vlines(0.0, 0.0, 1.1, 'k')
p.legend(loc='upper left')
p.grid()
p.draw()

p.savefig(paper1_figures_path("mask_edge_radial_profiles.pdf"))
p.savefig(paper1_figures_path("mask_edge_radial_profiles.png"))

# raw_input("Next plot?")
p.clf()

# Show the total number of elements in each distance bin

p.semilogy(bin_centers, binned_elements, 'bD-')
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Number of pixels")
p.grid()

p.savefig(paper1_figures_path("mask_edge_radial_profiles_numbin.pdf"))
p.savefig(paper1_figures_path("mask_edge_radial_profiles_numbin.png"))

# raw_input("Next plot?")
p.clf()

# Now investigate the significance of the distance correlations.
# Randomize the order of the CO and HI intensities.

hi_rand_vals = \
    binned_statistic(all_dists, np.random.permutation(all_vals_hi),
                     bins=bins,
                     statistic=np.nanmean)[0]

co_rand_vals = \
    binned_statistic(all_dists, np.random.permutation(all_vals_co),
                     bins=bins,
                     statistic=np.nanmean)[0]

p.plot(bin_centers, hi_rand_vals / np.nanmax(hi_rand_vals), 'bD-',
       label="HI")
p.plot(bin_centers, co_rand_vals / np.nanmax(co_rand_vals), 'ro-',
       label="CO(2-1)")
# p.xlim([0.0, 200])
p.ylim([0.0, 1.1])
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Normalized Intensity")
p.legend(loc='upper left')
p.grid()
p.draw()

p.savefig(paper1_figures_path("mask_edge_radial_profiles_randbin.pdf"))
p.savefig(paper1_figures_path("mask_edge_radial_profiles_randbin.png"))

# raw_input("Next plot?")
p.clf()

# Compare the CDFs of the intensities within the masks to demonstrate CO
# is not colocated with all of the HI

pos_hi = all_vals_hi[all_dists > 0]
pos_co = all_vals_co[all_dists > 0]

# Using argsort since
p.plot(np.sort(pos_hi), np.cumsum(np.sort(pos_hi)) / np.sum(pos_hi), "b-",
       label="HI")
p.plot(np.sort(pos_hi), np.cumsum(pos_co[np.argsort(pos_hi)]) / np.sum(pos_co),
       "g--", label="CO")
p.legend(loc='upper left')
p.grid()
p.ylim([-0.05, 1.05])
p.ylabel("CDF")
p.xlabel("HI Integrated Intensity (K km/s)")
p.draw()

p.savefig(paper1_figures_path("inmask_hi_co_cdfs.pdf"))
p.savefig(paper1_figures_path("inmask_hi_co_cdfs.png"))
p.close()

# Perform the same analysis split up into radial bins
dr = 500 * u.pc

max_radius = max_radius.to(u.pc)

nbins = np.int(np.floor(max_radius / dr))

inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)

Nrows = 4
Ncols = 3

p.figure(1, figsize=(12, 20)).clf()

fig, ax = p.subplots(Nrows, Ncols,
                     sharex=True,
                     sharey=True, num=1)

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)

fig.text(0.5, 0.02, 'Distance from mask edge (pc)', ha='center')
fig.text(0.04, 0.5, 'Normalized Intensity', va='center', rotation='vertical')

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    r, c = np.unravel_index(ctr, (Nrows, Ncols))

    idx = np.logical_and(all_radii >= r0.value,
                         all_radii < r1.value)

    hi_vals, bin_edges, bin_num = \
        binned_statistic(all_dists[idx], all_vals_hi[idx],
                         bins=bins,
                         statistic=np.nanmean)

    co_vals, bin_edges, bin_num = \
        binned_statistic(all_dists[idx], all_vals_co[idx],
                         bins=bins,
                         statistic=np.nanmean)

    binned_elements = \
        binned_statistic(all_dists[idx], np.ones_like(all_dists)[idx],
                         bins=bins,
                         statistic=np.sum)[0]

    hi_samps = np.zeros((niters, len(bin_centers)))
    co_samps = np.zeros((niters, len(bin_centers)))
    print("Bootstrapping")
    for i in ProgressBar(niters):
        hi_samps[i] = \
            binned_statistic(all_dists[idx],
                             np.random.permutation(all_vals_hi[idx]),
                             bins=bins,
                             statistic=np.nanmean)[0]

        co_samps[i] = \
            binned_statistic(all_dists[idx],
                             np.random.permutation(all_vals_co[idx]),
                             bins=bins,
                             statistic=np.nanmean)[0]

    # Take the stds in the distribution for each bin
    hi_errs = np.nanstd(hi_samps, axis=0)
    co_errs = np.nanstd(co_samps, axis=0)

    ax[r, c].errorbar(bin_centers, hi_vals / np.nanmax(hi_vals),
                      yerr=hi_errs / np.nanmax(hi_vals),
                      fmt="-",
                      color="b", label="HI")
    ax[r, c].errorbar(bin_centers, co_vals / np.nanmax(co_vals),
                      yerr=co_errs / np.nanmax(co_vals),
                      fmt="--",
                      color="r", label="CO(2-1)")
    ax[r, c].annotate("{0} to {1}".format(r0.to(u.kpc).value, r1.to(u.kpc)),
                      xy=(-360, 0.65),
                      color='k',
                      fontsize=14)

    ax[r, c].set_ylim([0.0, 1.1])
    # ax[r, c].set_xlabel("Distance from mask edge (pc)")
    # ax[r, c].set_ylabel("Normalized Intensity")
    # p.title("Radii {} to {}".format(r0, r1))
    ax[r, c].vlines(0.0, 0.0, 1.1, 'k')
    if ctr == 0:
        ax[r, c].legend(loc='upper left', fontsize=14)
    ax[r, c].grid()

    ax[r, c].set_xticklabels(ax[r, c].xaxis.get_majorticklabels(),
                             rotation=45)

fig.savefig(paper1_figures_path("mask_edge_radial_profiles_byradius.pdf"))
fig.savefig(paper1_figures_path("mask_edge_radial_profiles_byradius.png"))

p.close()
