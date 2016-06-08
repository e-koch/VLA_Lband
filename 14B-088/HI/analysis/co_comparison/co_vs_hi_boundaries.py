from basics import BubbleFinder2D
from basics.utils import sig_clip
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import astropy.units as u
import matplotlib.pyplot as p
import scipy.ndimage as nd
from scipy.stats import binned_statistic
import numpy as np
from skimage.morphology import medial_axis
from galaxies import Galaxy
from reproject import reproject_interp
from skimage.segmentation import find_boundaries
from astropy.utils.console import ProgressBar

'''
Calculate the intensities of HI and CO as a function of distance from the
edges of the adap. thresh mask.

To avoid a whole cube regrid, I use the rotation subtracted cubes and extract
zeroth moments from spectral slabs.

'''

# Plot a bunch
verbose = False

gal = Galaxy("M33")

# Load in the rotation subtracted cubes
hi_cube = SpectralCube.read("/home/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits")
co_cube = SpectralCube.read("/media/eric/Data_3/M33/co21/m33.co21_iram.rotsub.fits")

start_vel = - 30 * u.km / u.s
end_vel = 30 * u.km / u.s
del_vel = 3. * u.km / u.s

vels = np.arange(start_vel.value, end_vel.value + del_vel.value,
                 del_vel.value) * u.km / u.s

# Get the radius array so we can cut to where the CO data is valid
radii = gal.radius(header=hi_cube[0].header)
max_radius = 6.6 * u.kpc


all_dists = []
all_vals_hi = []
all_vals_co = []
edge_masks = []
skeletons = []

# Estimate the noise level in an equivalent slab
hi_mom0 = hi_cube.spectral_slab(-180 * u.km / u.s, -183 * u.km / u.s).moment0()
sigma = sig_clip(hi_mom0.value, nsig=10) * hi_cube.beam.jtok(1.414 * u.GHz).value

for i, (end, start) in enumerate(ProgressBar(zip(vels[1:], vels[:-1]))):

    hi_slab = hi_cube.spectral_slab(start, end)
    hi_mom0 = hi_slab.moment0() * hi_cube.beam.jtok(1.414 * u.GHz) / u.Jy

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
    radial_cut = radii <= 6.6 * u.kpc
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
    all_vals_hi.extend(list(hi_mom0.value[radial_cut] / 1000.))
    all_vals_co.extend(list(co_mom0.value[radial_cut] / 1000.))

    if verbose:
        print("Velocities: {} to {}".format(start, end))
        p.subplot(121)
        p.imshow(np.arctan(hi_mom0.value / np.nanpercentile(hi_mom0.value, 85)), origin='lower')
        # p.contour(skeleton, colors='b')
        p.contour(edge_mask, colors='g')

        p.subplot(122)
        p.imshow(hole_mask, origin='lower')
        # p.imshow(co_mom0.value, origin='lower')
        # p.contour(skeleton, colors='b')
        p.contour(edge_mask, colors='g')
        p.draw()

        raw_input("Next plot?")
        p.clf()

all_dists = np.array(all_dists)
all_vals_hi = np.array(all_vals_hi)
all_vals_co = np.array(all_vals_co)

# Now bin all of the distances against the HI and CO intensities.
bins = np.arange(-30, 30, 1)
hi_vals, bin_edges, bin_num = \
    binned_statistic(all_dists, all_vals_hi,
                     bins=bins,
                     statistic=np.nanmean)

hi_stds = \
    binned_statistic(all_dists, all_vals_hi,
                     bins=bins,
                     statistic=np.nanstd)[0]

co_vals, bin_edges, bin_num = \
    binned_statistic(all_dists, all_vals_co,
                     bins=bins,
                     statistic=np.nanmean)

co_stds = \
    binned_statistic(all_dists, all_vals_hi,
                     bins=bins,
                     statistic=np.nanstd)[0]

binned_elements = \
    binned_statistic(all_dists, np.ones_like(all_dists), bins=bins,
                     statistic=np.sum)[0]

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width / 2

# Convert the bin_centers to pc
pixscale = \
    hi_mom0.header['CDELT2'] * (np.pi / 180.) * 0.84e6  # * u.pc

bin_centers *= pixscale

# p.errorbar(bin_centers, hi_vals / np.nanmax(hi_vals),
#            yerr=hi_stds / np.nanmax(hi_vals), fmt="D",
#            color="b")
# p.errorbar(bin_centers, co_vals / np.nanmax(co_vals),
#            yerr=co_stds / np.nanmax(co_vals), fmt="o",
#            color="g")

p.plot(bin_centers, hi_vals / np.nanmax(hi_vals), 'bD-',
       label="HI")
p.plot(bin_centers, co_vals / np.nanmax(co_vals), 'ro-',
       label="CO(2-1)")
# p.xlim([0.0, 200])
p.ylim([0.0, 1.1])
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Normalized Intensity")
p.legend()
p.grid()
p.draw()

raw_input("Next plot?")
p.clf()

# Show the total number of elements in each distance bin

p.semilogy(bin_centers, binned_elements, 'bD-')
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Number of pixels")
p.grid()

raw_input("Next plot?")
p.clf()

# Now investigate the significance of the distance correlations.
# Randomize the order of the CO and HI intensities.

np.random.seed(34678953)

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
