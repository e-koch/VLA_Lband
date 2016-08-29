
import os
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.utils.console import ProgressBar
import numpy as np
from astropy.io import fits
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as p
from glob import glob

from basics import BubbleFinder, Bubble2D
from basics.utils import sig_clip

'''
Create the bubble catalogue of M33 with the 14B-088 map
'''

data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"
# data_path = "/lustre/home/ekoch/m33/"

# Prefix for save files
name = "M33_14B-088"

cube = \
    SpectralCube.read(os.path.join(data_path,
                                   "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits"))

lwidth = fits.getdata(os.path.join(data_path,
    "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.lwidth.fits"))

sigma = sig_clip(cube[-1].value, nsig=10)

# Create a cruddy linewidth map. This will need to be improved, but that can
# be done with just the saved bubble objects
# lwidth = cube.with_mask(cube > 2 * sigma * u.Jy).linewidth_sigma()
lwidth = lwidth * u.m / u.s

galaxy_props = {"center_coord":
                SkyCoord(23.461667, 30.660194,
                         unit=(u.deg, u.deg), frame='fk5'),
                "inclination": 56. * u.deg,
                "position_angle": 201. * u.deg,
                "scale_height": 100. * u.pc}


scales = 3. * np.arange(1, 15, np.sqrt(2))

# The first few channels have some emission in them. It's partially M33,
# partially galactic HI. So estimate the noise level from the last channel
# which is pretty much noise only. It gives 1.8 mJy/bm, which is right on.
bub_find = BubbleFinder(cube[:, 746:903, 556:747], keep_threshold_mask=True,
                        # empty_channel=cube.shape[0] - 1,  # Given above
                        sigma=sigma,
                        galaxy_props=galaxy_props,
                        distance=0.84 * u.Mpc)

# Save to separate folder for the sub-region
sep_folder = os.path.join(data_path, "sub_region_bubbles")
try:
    os.mkdir(sep_folder)
except OSError:
    pass

# Setup folder for saving the 2D regions
folder_2D = os.path.join(sep_folder, "bubbles_2D")
try:
    os.mkdir(folder_2D)
except OSError:
    pass

# Setup folder for the 3D bubbles and other output
# folder = os.path.expanduser("~/MyRAID/M33/14B-088/HI/full_imaging/bubbles/")
# folder = os.path.join(data_path, "bubbles/")
# try:
#     os.mkdir(folder)
# except OSError:
#     pass

# Load in existing 2D regions or not?
load_twods = False

# Requiring minimum of 15 channels, similar to the 3 channel requirement with
# LITTLE THINGS (15 * 0.2 km/s = 3 km/s; 3 * 1.2 km/s = 3.6 km/s)
if load_twods:
    # Load the mask
    mask = np.load(os.path.join(sep_folder, "{}_bubble_mask.npy".format(name)))

    # Load the 2D regions.
    files = glob(os.path.join(data_path, "bubbles_2D", "*.pkl"))

    twod_regions = []
    print("Loading in 2D regions.")
    for f in ProgressBar(files):
        twod_regions.append(Bubble2D.load_bubble(f))

    bub_find.get_bubbles(verbose=True, overlap_frac=0.5, multiprocess=True,
                         refit=False, nsig=1.5, min_corr=0.7, min_overlap=0.7,
                         global_corr=0.5, min_channels=15, nprocesses=None,
                         scales=scales, cube_linewidth=lwidth,
                         save_regions=False, twod_regions=twod_regions)
else:
    bub_find.get_bubbles(verbose=True, overlap_frac=0.5, multiprocess=True,
                         refit=False, nsig=1.5, min_corr=0.6, min_overlap=0.6,
                         global_corr=0.3, min_channels=15, nprocesses=None,
                         scales=scales, cube_linewidth=lwidth,
                         save_regions=True,
                         save_region_path=folder_2D)

# Save all of the bubbles
# bub_find.save_bubbles(folder=folder, name=name)

# catalog = bub_find.to_catalog()

# catalog.write_table(os.path.join(folder, "bubble_catalog.ecsv"))

# # Save the mask
# save_mask = False
# if save_mask:
#     np.save(os.path.join(folder, "{}_bubble_mask.npy".format(name)),
#             bub_find.mask)


# # Show the outlines of the bubbles
# fig = p.figure()
# ax = fig.add_subplot(111)
# ax = bub_find.visualize_bubbles(show=False, ax=ax,
#                                 plot_twoD_shapes=False)
# fig.savefig(os.path.join(folder, "{}_mom0_bubbles.pdf".format(name)))
# p.close()
# # With the 2D regions
# fig = p.figure()
# ax = fig.add_subplot(111)
# ax = bub_find.visualize_bubbles(show=False, ax=ax,
#                                 plot_twoD_shapes=True)
# fig.savefig(os.path.join(folder,
#                          "{}_mom0_bubbles_w_twoD.pdf".format(name)))
# p.close()

# # Individual channel maps
# output_folder_chans = os.path.join(folder, "channel_maps")
# try:
#     os.mkdir(output_folder_chans)
# except OSError:
#     pass
# bub_find.visualize_channel_maps(show_mask_contours=True, plot_unclustered=True,
#                                 save=True, save_path=output_folder_chans,
#                                 save_name=name)

# # Pruning off the arm regions. Basics cannot yet distinguish between holes
# # and spiral arm structure
# remove_bubbles = []

