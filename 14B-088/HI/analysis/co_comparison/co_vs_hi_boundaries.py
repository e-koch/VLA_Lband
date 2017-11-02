from basics import BubbleFinder2D
from basics.utils import sig_clip
from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as p
import scipy.ndimage as nd
from scipy.stats import binned_statistic
import numpy as np
from skimage.segmentation import find_boundaries
from skimage.morphology import medial_axis
from astropy.utils.console import ProgressBar
from corner import hist2d
import os
from os.path import join as osjoin
import seaborn as sb

from paths import (fourteenB_HI_data_wGBT_path,
                   iram_co21_14B088_data_path,
                   allfigs_path)
from constants import hi_freq
from galaxy_params import gal_feath as gal
from plotting_styles import (default_figure, twocolumn_figure,
                             onecolumn_figure,
                             twocolumn_twopanel_figure)


default_figure()

'''
Calculate the intensities of HI and CO as a function of distance from the
edges of the adap. thresh mask.

'''

fig_path = allfigs_path("co_vs_hi")
if not os.path.exists(fig_path):
    os.mkdir(fig_path)


np.random.seed(34678953)

# Plot a bunch
verbose = False
# slicer = (slice(825, 1033), slice(360, 692))
slicer = (slice(None), slice(None))

# Load in the rotation subtracted cubes
hi_cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_to_co/M33_14B-088_HI.clean.image.GBT_feathered.2.6kms.fits"))
co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.fits"))

co_cube = co_cube.spectral_slab(hi_cube.spectral_extrema[0],
                                hi_cube.spectral_extrema[1])

# Skip the first 7 channels and the last 15 channels
vels = co_cube.spectral_axis.to(u.km / u.s)[7:-15]

# Get the radius array so we can cut to where the CO data is valid
radii = gal.radius(header=hi_cube[0].header)
max_radius = 6.0 * u.kpc

all_dists = []
all_radii = []
all_vals_hi = []
all_vals_co = []
edge_masks = []
skeleton_dists = []
skeleton_dists_pix = []
skeleton_widths = []

masks = []

hi_beam = hi_cube.beam

# Estimate the noise level in an equivalent slab
hi_mom0 = hi_cube[-1]
sigma = sig_clip(hi_mom0.value, nsig=10) * \
    hi_beam.jtok(hi_freq).value

# Skip the first 7 channels
i_offset = 7

for i, vel in enumerate(ProgressBar(vels)):

    hi_chan = hi_cube[i + i_offset] * hi_beam.jtok(hi_freq) / u.Jy
    co_chan = co_cube[i + i_offset]

    # Need a mask from the HI
    # Adjust the sigma in a single channel to the moment0 in the slab
    # sigma = 0.00152659 * hi_slab.shape[0] * \
    #     np.abs((hi_slab.spectral_axis[1] - hi_slab.spectral_axis[0]).value)

    bub = BubbleFinder2D(hi_chan, auto_cut=False, sigma=sigma)
    bub.create_mask(bkg_nsig=5, region_min_nsig=10, mask_clear_border=False)

    skeleton, dists = medial_axis(~bub.mask, return_distance=True)
    skeleton_dists.append(skeleton * dists)

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
    all_vals_hi.extend(list(hi_chan.value[radial_cut]))
    all_vals_co.extend(list(co_chan.value[radial_cut]))

    # Track the width of the mask
    skeleton_widths.extend(list((skeleton * dists)[radial_cut]))
    # Also record all of the distances from the centre of the skeletons
    skeleton_dists_pix.extend(list(nd.distance_transform_edt(~skeleton)[radial_cut]))

    masks.append(~bub.mask)

    if verbose:
        print("Velocity: {}".format(vel))

        fig, ax = p.subplots(1, 2, sharex=True, sharey=True,
                             subplot_kw=dict(projection=hi_chan.wcs))

        ax[0].imshow(np.arctan(hi_chan[slicer].value /
                               np.nanpercentile(hi_chan.value, 85)),
                     origin='lower', vmin=0)
        # p.contour(skeleton, colors='b')
        ax[0].contour(edge_mask[slicer], colors='g')

        # p.imshow(hole_mask, origin='lower')
        ax[1].imshow(co_chan.value[slicer],
                     origin='lower', vmin=0.01, vmax=0.1)
        # p.imshow(co_chan.value, origin='lower')
        # p.contour(skeleton, colors='b')
        ax[1].contour(edge_mask[slicer], colors='g')
        lat = ax[1].coords[1]
        lat.set_ticklabel_visible(False)

        p.draw()
        raw_input("Next plot?")
        p.clf()

all_dists = np.array(all_dists)
all_radii = np.array(all_radii)
all_vals_hi = np.array(all_vals_hi)
all_vals_co = np.array(all_vals_co)

skeleton_widths = np.array(skeleton_widths)
skeleton_dists_pix = np.array(skeleton_dists_pix)

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

onecolumn_figure()

p.errorbar(bin_centers, hi_vals / np.nanmax(hi_vals),
           yerr=hi_errs / np.nanmax(hi_vals), fmt="D-",
           label="HI")
p.errorbar(bin_centers, co_vals / np.nanmax(co_vals),
           yerr=co_errs / np.nanmax(co_vals), fmt="o-",
           label="CO(2-1)")

# p.plot(bin_centers, hi_vals / np.nanmax(hi_vals), 'bD-',
#        label="HI")
# p.plot(bin_centers, co_vals / np.nanmax(co_vals), 'ro-',
#        label="CO(2-1)")
# p.xlim([0.0, 200])
p.ylim([0.0, 1.1])
p.xlim([-400, 220])
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Normalized Intensity")
p.vlines(0.0, 0.0, 1.1, 'k')
p.legend(loc='upper left', frameon=True)
p.grid()
p.tight_layout()
p.draw()

p.savefig(osjoin(fig_path, "mask_edge_radial_profiles.pdf"))
p.savefig(osjoin(fig_path, "mask_edge_radial_profiles.png"))

# raw_input("Next plot?")
p.close()

# Show the total number of elements in each distance bin

p.semilogy(bin_centers, binned_elements, 'D-')
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Number of pixels")
p.grid()
p.tight_layout()

p.savefig(osjoin(fig_path, "mask_edge_radial_profiles_numbin.pdf"))
p.savefig(osjoin(fig_path, "mask_edge_radial_profiles_numbin.png"))

# raw_input("Next plot?")
p.close()

# Now investigate the significance of the distance correlations.
# Randomize the order of the CO and HI intensities.

# hi_rand_vals = \
#     binned_statistic(all_dists, np.random.permutation(all_vals_hi),
#                      bins=bins,
#                      statistic=np.nanmean)[0]

# co_rand_vals = \
#     binned_statistic(all_dists, np.random.permutation(all_vals_co),
#                      bins=bins,
#                      statistic=np.nanmean)[0]

# p.plot(bin_centers, hi_rand_vals / np.nanmax(hi_rand_vals), 'bD-',
#        label="HI")
# p.plot(bin_centers, co_rand_vals / np.nanmax(co_rand_vals), 'ro-',
#        label="CO(2-1)")
# # p.xlim([0.0, 200])
# p.ylim([0.0, 1.1])
# p.xlabel("Distance from mask edge (pc)")
# p.ylabel("Normalized Intensity")
# p.legend(loc='upper left')
# p.grid()
# p.draw()

# p.savefig(paper1_figures_path("mask_edge_radial_profiles_randbin.pdf"))
# p.savefig(paper1_figures_path("mask_edge_radial_profiles_randbin.png"))

# # raw_input("Next plot?")
# p.close()

# Compare the CDFs of the intensities within the masks to demonstrate CO
# is not colocated with all of the HI

pos_hi = all_vals_hi[all_dists > 0]
pos_co = all_vals_co[all_dists > 0]

onecolumn_figure()

p.plot(np.sort(pos_hi), np.cumsum(np.sort(pos_hi)) / np.sum(pos_hi), "-",
       label="HI")
p.plot(np.sort(pos_hi), np.cumsum(pos_co[np.argsort(pos_hi)]) / np.sum(pos_co),
       "--", label="CO")
p.legend(loc='upper left', frameon=True)
p.grid()
p.ylim([-0.05, 1.05])
p.ylabel("CDF")
p.xlabel("HI Intensity (K)")
p.tight_layout()

p.savefig(osjoin(fig_path, "inmask_hi_co_cdfs.pdf"))
p.savefig(osjoin(fig_path, "inmask_hi_co_cdfs.png"))
p.close()

# Perform the same analysis split up into radial bins
dr = 500 * u.pc

max_radius = max_radius.to(u.pc)

nbins = np.int(np.floor(max_radius / dr))

inneredge = np.linspace(0, max_radius - dr, nbins)
outeredge = np.linspace(dr, max_radius, nbins)

Nrows = 4
Ncols = 3

twocolumn_figure()
p.figure(1, figsize=(8.4, 11)).clf()

fig, ax = p.subplots(Nrows, Ncols,
                     sharex=True,
                     sharey=True, num=1)

fig.text(0.5, 0.04, 'Distance from mask edge (pc)', ha='center')
fig.text(0.04, 0.5, 'Normalized Intensity', va='center', rotation='vertical')

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)

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
                      label="HI")
    ax[r, c].errorbar(bin_centers, co_vals / np.nanmax(co_vals),
                      yerr=co_errs / np.nanmax(co_vals),
                      fmt="--",
                      label="CO(2-1)")
    ax[r, c].annotate("{0} to {1}".format(r0.to(u.kpc).value, r1.to(u.kpc)),
                      xy=(-360, 0.65),
                      color='k',
                      fontsize=12,
                      bbox={"boxstyle": "square", "facecolor": "w"})

    ax[r, c].set_ylim([0.0, 1.1])
    ax[r, c].set_xlim([-400, 220])
    # ax[r, c].set_xlabel("Distance from mask edge (pc)")
    # ax[r, c].set_ylabel("Normalized Intensity")
    # p.title("Radii {} to {}".format(r0, r1))
    ax[r, c].vlines(0.0, 0.0, 1.1, 'k')
    if ctr == 0:
        ax[r, c].legend(loc='upper left', fontsize=12, frameon=True)
    ax[r, c].grid()

# for r in range(Nrows):
#     for c in range(Ncols):
#         if r == Nrows - 1:
#             ax[r, c].set_xticklabels(ax[r, c].xaxis.get_majorticklabels(),
#                                      rotation=45)

fig.savefig(osjoin(fig_path, "mask_edge_radial_profiles_byradius.pdf"))
fig.savefig(osjoin(fig_path, "mask_edge_radial_profiles_byradius.png"))

p.close()

onecolumn_figure()

# Is the variation being driven by a change in the width of the regions?
bins = np.arange(0, 6.5, 0.5) * 1000
dists, bin_edges, bin_num = \
    binned_statistic(all_radii[skeleton_widths > 0],
                     skeleton_widths[skeleton_widths > 0],
                     bins=bins, statistic=np.mean)
dist_std = \
    binned_statistic(all_radii[skeleton_widths > 0],
                     skeleton_widths[skeleton_widths > 0],
                     bins=bins, statistic=np.std)[0]

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width / 2

ang_conv = (hi_mom0.header["CDELT2"] * u.deg).to(u.arcsec)
phys_conv = ang_conv.to(u.rad).value * 840e3 * u.pc

p.plot(all_radii[skeleton_widths > 0] / 1000.,
       skeleton_widths[skeleton_widths > 0] * phys_conv.value,
       'ko', alpha=0.1, ms=3.0, zorder=-1)

p.errorbar(bin_centers / 1000., dists * phys_conv.value,
           yerr=dist_std * phys_conv.value,
           marker='', linestyle='-', linewidth=2, elinewidth=2)
p.ylabel("Width of mask regions (pc)")
p.xlabel("Radius (kpc)")

p.tight_layout()

p.savefig(osjoin(fig_path, "mask_width_byradius.pdf"))
p.savefig(osjoin(fig_path, "mask_width_byradius.png"))
p.close()

# Let's take some other views of this data, while we're at it.

twocolumn_twopanel_figure()

fig, ax = p.subplots(1, 2, sharex=True)

hist2d(skeleton_widths[skeleton_widths > 0] * phys_conv.value,
       all_vals_co[skeleton_widths > 0], bins=10,
       ax=ax[1], data_kwargs={"alpha": 0.6})
ax[1].set_xlabel("Mask Width (pc)")
ax[1].set_ylabel(r"CO Intensity (K)")
ax[1].grid()

hist2d(skeleton_widths[skeleton_widths > 0] * phys_conv.value,
       all_vals_hi[skeleton_widths > 0], bins=10,
       ax=ax[0], data_kwargs={"alpha": 0.6})
ax[0].set_xlabel("Mask Width (pc)")
ax[0].set_ylabel(r"HI Intensity (K)")
ax[0].grid()

p.tight_layout()

p.savefig(osjoin(fig_path, "mask_widthvsintensity.pdf"))
p.savefig(osjoin(fig_path, "mask_widthvsintensity.png"))
p.close()

# HI vs. CO with all skeleton distances (not just on the skeleton like above)
bins = np.arange(0, 31, 1)

selector_pts = np.logical_and(all_vals_co > 0,
                              np.logical_and(all_vals_hi > 0,
                                             skeleton_dists_pix < 30))

hi_mean, bin_edges, bin_num = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_hi[selector_pts],
                     bins=bins, statistic=np.mean)
hi_85 = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_hi[selector_pts],
                     bins=bins, statistic=lambda x: np.percentile(x, 85))[0]
hi_15 = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_hi[selector_pts],
                     bins=bins, statistic=lambda x: np.percentile(x, 15))[0]

co_mean = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_co[selector_pts],
                     bins=bins, statistic=np.mean)[0]
co_85 = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_co[selector_pts],
                     bins=bins, statistic=lambda x: np.percentile(x, 85))[0]
co_15 = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_co[selector_pts],
                     bins=bins, statistic=lambda x: np.percentile(x, 15))[0]

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width / 2

twocolumn_twopanel_figure()

fig, ax = p.subplots(1, 2, sharex=True)

ax[1].plot(skeleton_dists_pix[selector_pts] * phys_conv.value,
           all_vals_co[selector_pts], 'ko', ms=2.0, alpha=0.6,
           rasterized=True, zorder=-1)
ax[1].set_xlabel("Distance from Mask Centre (pc)")
ax[1].set_ylabel(r"CO Intensity (K)")
ax[1].grid()

ax[0].plot(skeleton_dists_pix[selector_pts] * phys_conv.value,
           all_vals_hi[selector_pts], 'ko', ms=2.0, alpha=0.6,
           rasterized=True, zorder=-1)
ax[0].set_xlabel("Distance from Mask Centre (pc)")
ax[0].set_ylabel(r"HI Intensity (K)")
ax[0].grid()

p.tight_layout()

p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist.pdf"))
p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist.png"))
p.close()

onecolumn_figure()

col_pal = sb.color_palette()

ax = p.subplot(111)
pl1 = ax.plot(bin_centers * phys_conv.value, hi_mean, "D-",
              label='HI')
# ax.errorbar(bin_centers * phys_conv.value, hi_mean,
#             yerr=[hi_mean - hi_15, hi_85 - hi_mean], color='b',
#             marker='D', linestyle='-', linewidth=2, elinewidth=2)
ax.set_xlabel("Distance from Mask Centre (pc)")
ax.set_ylabel(r"HI Mean Intensity (K)")

ax_2 = ax.twinx()
pl2 = ax_2.plot(bin_centers * phys_conv.value, co_mean, "o--",
                label='CO', color=col_pal[1])
# ax_2.errorbar(bin_centers * phys_conv.value, co_mean,
#               yerr=[co_mean - co_15, co_85 - co_mean], color='r',
#               marker='o', linestyle='--', linewidth=2, elinewidth=2)
ax_2.set_ylabel(r"CO Mean Intensity (K)")
# align_yaxis(ax, 0, ax_2, 0)

pls = pl1 + pl2
labs = [l.get_label() for l in pls]
ax.legend(pls, labs, frameon=True)

ax.grid()

p.tight_layout()

p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_mean.pdf"))
p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_mean.png"))
p.close()


default_figure()
