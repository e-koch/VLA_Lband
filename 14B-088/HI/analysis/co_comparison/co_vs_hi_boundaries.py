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
from aplpy import FITSFigure
from astropy.io import fits
from astropy.table import Table, Column

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

# Make a figure from one of the channels to highlight the mask shape
twocolumn_twopanel_figure()

mpl_fig = p.figure()

spatial_slice = (slice(720, 900), slice(350, 750))

fig = FITSFigure((hi_cube[47][spatial_slice] * hi_beam.jtok(hi_freq).value).hdu,
                 figure=mpl_fig)
fig.show_grayscale(invert=True, vmin=None, vmax=80, stretch='sqrt')
fig.add_colorbar()
fig.colorbar.set_axis_label_text("HI Intensity (K)")
fig.show_contour(co_cube[47][spatial_slice].hdu, cmap='autumn',
                 levels=[0.05, 0.1, 0.2, 0.3])
fig.show_contour(fits.PrimaryHDU(masks[40].astype(int), hi_cube[0].header),
                 colors=[sb.color_palette()[-1]], levels=[0.5])
fig.hide_axis_labels()

p.tight_layout()

fig.savefig(osjoin(fig_path, "mask_edge_img_vel_minus196.pdf"))
fig.savefig(osjoin(fig_path, "mask_edge_img_vel_minus196.png"))
fig.close()

# Make a version with the skeleton instead of the mask
mpl_fig = p.figure()

spatial_slice = (slice(720, 900), slice(350, 750))

fig = FITSFigure((hi_cube[47][spatial_slice] * hi_beam.jtok(hi_freq).value).hdu,
                 figure=mpl_fig)
fig.show_grayscale(invert=True, vmin=None, vmax=80, stretch='sqrt')
fig.add_colorbar()
fig.colorbar.set_axis_label_text("HI Intensity (K)")
fig.show_contour(fits.PrimaryHDU(medial_axis(masks[40]).astype(int),
                                 hi_cube[0].header),
                 colors=[sb.color_palette()[-1]], levels=[0.5])
fig.show_contour(co_cube[47][spatial_slice].hdu, cmap='autumn',
                 levels=[0.05, 0.1, 0.2, 0.3])
fig.hide_axis_labels()

p.tight_layout()

fig.savefig(osjoin(fig_path, "mask_skel_img_vel_minus196.pdf"))
fig.savefig(osjoin(fig_path, "mask_skel_img_vel_minus196.png"))
fig.close()

# Now bin all of the distances against the HI and CO intensities.
bins = np.arange(-30, 30, 1)
hi_vals, bin_edges, bin_num = \
    binned_statistic(all_dists, all_vals_hi,
                     bins=bins,
                     statistic=np.mean)

co_vals, bin_edges, bin_num = \
    binned_statistic(all_dists, all_vals_co,
                     bins=bins,
                     statistic=np.mean)

binned_elements = \
    binned_statistic(all_dists, np.ones_like(all_dists), bins=bins,
                     statistic=np.sum)[0]

# Require that there be 100 points in each bin
bin_cutoff = binned_elements >= 100

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width / 2

bin_centers = bin_centers[bin_cutoff]
co_vals = co_vals[bin_cutoff]
hi_vals = hi_vals[bin_cutoff]

# Let's bootstrap to get errors in the distance bins
niters = 100
hi_samps = np.zeros((niters, len(bin_centers)))
co_samps = np.zeros((niters, len(bin_centers)))
print("Bootstrapping")
for i in ProgressBar(niters):
    hi_samps[i] = \
        binned_statistic(all_dists, np.random.permutation(all_vals_hi),
                         bins=bins,
                         statistic=np.mean)[0][bin_cutoff]

    co_samps[i] = \
        binned_statistic(all_dists, np.random.permutation(all_vals_co),
                         bins=bins,
                         statistic=np.mean)[0][bin_cutoff]

# Take the stds in the distribution for each bin
hi_errs = np.nanstd(hi_samps, axis=0)
co_errs = np.nanstd(co_samps, axis=0)

# Convert the bin_centers to pc
pixscale = \
    hi_mom0.header['CDELT2'] * (np.pi / 180.) * gal.distance.to(u.pc).value

bin_centers *= pixscale

onecolumn_figure()
cpal = sb.color_palette()

ax = p.subplot(111)

pl1 = ax.errorbar(bin_centers, hi_vals,
                  yerr=hi_errs, fmt="D-", color=cpal[0],
                  label="HI", drawstyle='steps-mid')
ax_2 = ax.twinx()
pl2 = ax_2.errorbar(bin_centers, co_vals * 1000.,
                    yerr=co_errs, fmt="o--",
                    color=cpal[1],
                    label="CO(2-1)", drawstyle='steps-mid')

pls = [pl1[0], pl2[0]]
labs = ["HI", "CO(2-1)"]
ax.legend(pls, labs, frameon=True)

ax.set_xlabel("Distance from Mask Edge (pc)")
ax.set_ylabel(r"HI Mean Intensity (K)")
ax_2.set_ylabel(r"CO Mean Intensity (mK)")

ax.axvline(0.0, color='k')
ax.grid()

p.tight_layout()

p.savefig(osjoin(fig_path, "mask_edge_radial_profiles.pdf"))
p.savefig(osjoin(fig_path, "mask_edge_radial_profiles.png"))

p.close()

# Show the total number of elements in each distance bin

p.semilogy(pixscale * (bin_edges[1:] - bin_width / 2),
           binned_elements, 'D-')
p.axhline(100, linestyle='--', color=sb.color_palette()[1])
p.xlabel("Distance from mask edge (pc)")
p.ylabel("Number of pixels")
p.grid()
p.tight_layout()

p.savefig(osjoin(fig_path, "mask_edge_radial_profiles_numbin.pdf"))
p.savefig(osjoin(fig_path, "mask_edge_radial_profiles_numbin.png"))

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

hi_vals_rad = []
hi_errs_rad = []
co_vals_rad = []
co_errs_rad = []

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    r, c = np.unravel_index(ctr, (Nrows, Ncols))

    idx = np.logical_and(all_radii >= r0.value,
                         all_radii < r1.value)

    hi_vals_bin, bin_edges, bin_num = \
        binned_statistic(all_dists[idx], all_vals_hi[idx],
                         bins=bins,
                         statistic=np.mean)

    co_vals_bin, bin_edges, bin_num = \
        binned_statistic(all_dists[idx], all_vals_co[idx],
                         bins=bins,
                         statistic=np.mean)

    binned_elements = \
        binned_statistic(all_dists[idx], np.ones_like(all_dists)[idx],
                         bins=bins,
                         statistic=np.sum)[0]

    bin_cutoff = binned_elements >= 30

    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width / 2

    bin_centers = bin_centers[bin_cutoff]
    co_vals_bin = co_vals_bin[bin_cutoff]
    hi_vals_bin = hi_vals_bin[bin_cutoff]

    hi_samps = np.zeros((niters, len(bin_centers)))
    co_samps = np.zeros((niters, len(bin_centers)))
    print("Bootstrapping")
    for i in ProgressBar(niters):
        hi_samps[i] = \
            binned_statistic(all_dists[idx],
                             np.random.permutation(all_vals_hi[idx]),
                             bins=bins,
                             statistic=np.mean)[0][bin_cutoff]

        co_samps[i] = \
            binned_statistic(all_dists[idx],
                             np.random.permutation(all_vals_co[idx]),
                             bins=bins,
                             statistic=np.mean)[0][bin_cutoff]

    # Take the stds in the distribution for each bin
    hi_errs_bin = np.nanstd(hi_samps, axis=0)
    co_errs_bin = np.nanstd(co_samps, axis=0)

    hi_vals_rad.append(hi_vals)
    hi_errs_rad.append(hi_errs)
    co_vals_rad.append(co_vals)
    co_errs_rad.append(co_errs)

    ax[r, c].errorbar(bin_centers * pixscale,
                      hi_vals_bin / np.nanmax(hi_vals_bin),
                      yerr=hi_errs_bin / np.nanmax(hi_vals_bin),
                      fmt="D-", drawstyle='steps-mid',
                      label="HI")
    ax[r, c].errorbar(bin_centers * pixscale,
                      co_vals_bin / np.nanmax(co_vals_bin),
                      yerr=co_errs_bin / np.nanmax(co_vals_bin),
                      fmt="o--", drawstyle='steps-mid',
                      label="CO(2-1)")
    ax[r, c].annotate("{0} to {1}".format(r0.to(u.kpc).value, r1.to(u.kpc)),
                      xy=(-360, -0.3),
                      color='k',
                      fontsize=11,
                      bbox={"boxstyle": "square", "facecolor": "w"})

    # ax[r, c].set_ylim([0.0, 1.1])
    # ax[r, c].set_xlim([-400, 220])
    # ax[r, c].set_xlabel("Distance from mask edge (pc)")
    # ax[r, c].set_ylabel("Normalized Intensity")
    # p.title("Radii {} to {}".format(r0, r1))
    ax[r, c].axvline(0.0, color='k')
    if ctr == 0:
        ax[r, c].legend(loc='upper left', fontsize=11, frameon=True)
    ax[r, c].grid()

ax[r, c].set_ylim([-0.4, 1.1])

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
bins = np.arange(0, 36, 1)

selector_pts = np.logical_and(all_vals_co > 0,
                              np.logical_and(all_vals_hi > 0,
                                             skeleton_dists_pix < 35))

hi_mean, bin_edges, bin_num = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_hi[selector_pts],
                     bins=bins, statistic=np.mean)

hi_std = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_hi[selector_pts],
                     bins=bins, statistic=np.std)[0]

co_mean = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_co[selector_pts],
                     bins=bins, statistic=np.mean)[0]

co_std = \
    binned_statistic(skeleton_dists_pix[selector_pts],
                     all_vals_co[selector_pts],
                     bins=bins, statistic=np.std)[0]

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width / 2

num_in_bins = np.bincount(bin_num)[1:]
# Num. indep't points divided by number of pixels in one beam.
num_indept = num_in_bins / 41.

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


# Find the HWHM points for the radial profiles
def find_perc_width(x, y, level=0.5):
    '''
    Return the equivalent Gaussian sigma based on the HWHM positions.
    '''
    from scipy.interpolate import InterpolatedUnivariateSpline

    assert (level > 0.) & (level < 1.)

    # Assume that the profile peaks at the centre and monotonically
    # decreases. This is true for our comparisons here.
    peak = y.max()
    bkg = np.mean(y[-5:])

    halfmax = (peak - bkg) * level

    # Model the spectrum with a spline
    # x values must be increasing for the spline, so flip if needed.
    interp1 = InterpolatedUnivariateSpline(x, y - halfmax - bkg, k=3)

    hwhm_points = interp1.roots()
    if len(hwhm_points) < 1:
        raise ValueError("Didn't find HWHM!")
    elif len(hwhm_points) > 1:
        hwhm_points = [min(hwhm_points)]

    return hwhm_points[0]


co_hwhm = find_perc_width(bin_centers, co_mean) * phys_conv.value
hi_hwhm = find_perc_width(bin_centers, hi_mean) * phys_conv.value

# Widths at the 25 and 75th percentiles between peak and bkg.
co_lowq = find_perc_width(bin_centers, co_mean, level=0.25) * phys_conv.value
hi_lowq = find_perc_width(bin_centers, hi_mean, level=0.25) * phys_conv.value
co_highq = find_perc_width(bin_centers, co_mean, level=0.75) * phys_conv.value
hi_highq = find_perc_width(bin_centers, hi_mean, level=0.75) * phys_conv.value


print(co_hwhm, hi_hwhm)
# (67.480769229465864, 100.43889752239386)
print(co_lowq, hi_lowq)
# (103.06346550608011, 169.52993127884227)
print(co_highq, hi_highq)
# (40.128517012921016, 56.468698633140313)

onecolumn_figure()

col_pal = sb.color_palette()

ax = p.subplot(111)
ax.axvline(hi_hwhm, linestyle='-', color=col_pal[0],
           linewidth=3, alpha=0.6)
ax.axvline(co_hwhm, linestyle='--', color=col_pal[1],
           linewidth=3, alpha=0.6)

pl1 = ax.errorbar(bin_centers * phys_conv.value, hi_mean,
                  yerr=hi_std / np.sqrt(num_indept), color=col_pal[0],
                  marker='D', linestyle='-', linewidth=2, elinewidth=2,
                  label='HI')
ax.set_xlabel("Distance from Skeleton (pc)")
ax.set_ylabel(r"HI Mean Intensity (K)")

ax_2 = ax.twinx()
pl2 = ax_2.errorbar(bin_centers * phys_conv.value, co_mean,
                    yerr=co_std / np.sqrt(num_indept), color=col_pal[1],
                    marker='o', linestyle='--', linewidth=2, elinewidth=2,
                    label='CO')
ax_2.set_ylabel(r"CO Mean Intensity (K)")

pls = [pl1[0], pl2[0]]
labs = ['HI', 'CO']
ax.legend(pls, labs, frameon=True)

ax.grid()

p.tight_layout()

p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_mean.pdf"))
p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_mean.png"))
p.close()

# Split the profiles by radius
hi_hwhms = []
co_hwhms = []
hi_lowqs = []
co_lowqs = []
hi_highqs = []
co_highqs = []

twocolumn_figure()
p.figure(1, figsize=(8.4, 11)).clf()

fig, ax = p.subplots(Nrows, Ncols,
                     sharex=True,
                     sharey=True, num=1)

fig.text(0.5, 0.04, "Distance from Skeleton (pc)", ha='center')
fig.text(0.04, 0.5, r"HI Mean Intensity (K)", va='center', rotation='vertical')
fig.text(0.96, 0.5, r"CO Mean Intensity (K)", va='center', rotation='vertical')

p.subplots_adjust(hspace=0.1,
                  wspace=0.1)

for ctr, (r0, r1) in enumerate(zip(inneredge,
                                   outeredge)):

    r, c = np.unravel_index(ctr, (Nrows, Ncols))

    idx = np.logical_and(all_radii >= r0.value,
                         all_radii < r1.value)

    idx = np.logical_and(selector_pts, idx)

    hi_mean_rad, bin_edges, bin_num = \
        binned_statistic(skeleton_dists_pix[idx],
                         all_vals_hi[idx],
                         bins=bins, statistic=np.mean)
    hi_std_rad = \
        binned_statistic(skeleton_dists_pix[idx],
                         all_vals_hi[idx],
                         bins=bins, statistic=np.std)[0]

    co_mean_rad = \
        binned_statistic(skeleton_dists_pix[idx],
                         all_vals_co[idx],
                         bins=bins, statistic=np.mean)[0]

    co_std_rad = \
        binned_statistic(skeleton_dists_pix[idx],
                         all_vals_co[idx],
                         bins=bins, statistic=np.std)[0]

    num_in_bins_rad = np.bincount(bin_num)[1:]
    # Num. indep't points divided by number of pixels in one beam.
    num_indept_rad = num_in_bins / 41.

    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width / 2

    co_hwhm_rad = find_perc_width(bin_centers, co_mean_rad) * phys_conv.value
    hi_hwhm_rad = find_perc_width(bin_centers, hi_mean_rad) * phys_conv.value

    co_lowq_rad = find_perc_width(bin_centers, co_mean_rad,
                                  level=0.25) * phys_conv.value
    hi_lowq_rad = find_perc_width(bin_centers, hi_mean_rad,
                                  level=0.25) * phys_conv.value
    co_highq_rad = find_perc_width(bin_centers, co_mean_rad,
                                   level=0.75) * phys_conv.value
    hi_highq_rad = find_perc_width(bin_centers, hi_mean_rad,
                                   level=0.75) * phys_conv.value

    co_hwhms.append(co_hwhm_rad)
    hi_hwhms.append(hi_hwhm_rad)
    co_lowqs.append(co_lowq_rad)
    hi_lowqs.append(hi_lowq_rad)
    co_highqs.append(co_highq_rad)
    hi_highqs.append(hi_highq_rad)

    ax[r, c].axvline(hi_hwhm_rad, linestyle='-', color=col_pal[0],
                     linewidth=3, alpha=0.6)
    ax[r, c].axvline(co_hwhm_rad, linestyle='--', color=col_pal[1],
                     linewidth=3, alpha=0.6)

    pl1 = ax[r, c].errorbar(bin_centers * pixscale,
                            hi_mean_rad,
                            yerr=hi_std_rad / np.sqrt(num_indept_rad),
                            fmt="D-", drawstyle='steps-mid',
                            label="HI")
    ax2 = ax[r, c].twinx()
    if c == 0:
        ax_twin = ax2
    else:
        ax2.get_shared_y_axes().join(ax2, ax_twin)

    if c < Ncols - 1:
        ax2.tick_params(labelright='off')

    pl2 = ax2.errorbar(bin_centers * pixscale,
                       co_mean_rad,
                       yerr=co_std_rad / np.sqrt(num_indept_rad),
                       fmt="o--", drawstyle='steps-mid',
                       label="CO(2-1)", color=col_pal[1])
    ax[r, c].annotate("{0} to {1}".format(r0.to(u.kpc).value, r1.to(u.kpc)),
                      xy=(170, 25),
                      color='k',
                      fontsize=11,
                      bbox={"boxstyle": "square", "facecolor": "w"})

    ax[r, c].axvline(0.0, color='k')
    if ctr == 0:
        pls = [pl1[0], pl2[0]]
        labs = ["HI", "CO"]
        ax[r, c].legend(pls, labs, loc='center right', fontsize=11,
                        frameon=True)
    ax[r, c].grid()

fig.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_mean_byradius.pdf"))
fig.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_mean_byradius.png"))

p.close()


# Look at the HWHM with radius

onecolumn_figure()
p.plot(inneredge.value / 1000. + 0.25, hi_lowqs, label='HI 25\%',
       drawstyle='steps-mid', color=col_pal[0], linestyle='--')
p.plot(inneredge.value / 1000. + 0.25, co_lowqs, label='CO 25\%',
       drawstyle='steps-mid', color=col_pal[1], linestyle='--')

p.plot(inneredge.value / 1000. + 0.25, hi_hwhms, label='HI 50\%',
       drawstyle='steps-mid', color=col_pal[0], linewidth=3)
p.plot(inneredge.value / 1000. + 0.25, co_hwhms, label='CO 50\%',
       drawstyle='steps-mid', color=col_pal[1], linewidth=3)

p.plot(inneredge.value / 1000. + 0.25, hi_highqs, label='HI 75\%',
       drawstyle='steps-mid', color=col_pal[0], linestyle='-.')
p.plot(inneredge.value / 1000. + 0.25, co_highqs, label='CO 75\%',
       drawstyle='steps-mid', color=col_pal[1], linestyle='-.')
p.grid()
p.xlim([-0.2, 8.5])
p.legend(frameon=True, loc='center right')
p.ylabel("Width (pc)")
p.xlabel("Galactocentric Radius (kpc)")
p.tight_layout()

p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_hwhm.pdf"))
p.savefig(osjoin(fig_path, "mask_intensity_vs_skeldist_hwhm.png"))
p.close()


# Save the table of HWHMs
radbin_centers = Column((inneredge.value + 250) / 1000., unit=u.kpc,
                        name='bin_cent')
hi_hwhms = Column(hi_hwhms, unit=u.pc, name='hi_hwhm')
co_hwhms = Column(co_hwhms, unit=u.pc, name='co_hwhm')
hi_lowqs = Column(hi_lowqs, unit=u.pc, name='hi_25')
co_lowqs = Column(co_lowqs, unit=u.pc, name='co_25')
hi_highqs = Column(hi_highqs, unit=u.pc, name='hi_75')
co_highqs = Column(co_highqs, unit=u.pc, name='co_75')

tab = Table([radbin_centers, hi_hwhms, co_hwhms, hi_lowqs, co_lowqs, hi_highqs, co_highqs])

tab.write(fourteenB_HI_data_wGBT_path("tables/skeleton_profile_radial_hwhm.fits", no_check=True))

default_figure()
