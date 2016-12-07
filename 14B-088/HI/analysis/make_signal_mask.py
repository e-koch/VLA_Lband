
from astropy.io import fits
import astropy.units as u
from spectral_cube import SpectralCube
from signal_id import RadioMask, Noise
from scipy import ndimage as nd
from scipy.signal import medfilt
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.utils.console import ProgressBar
from astropy.stats import mad_std
import skimage.morphology as mo
import numpy as np
import skimage.morphology as mo
from spectral_cube.cube_utils import average_beams
from itertools import groupby, chain
from operator import itemgetter
import matplotlib.pyplot as p

from basics.utils import sig_clip

from analysis.paths import fourteenB_HI_data_path
from analysis.constants import cube_name, mask_name

'''
Create a signal mask for the 14B-088 cube
'''


def round_up_to_odd(f):
    return np.ceil(f) // 2 * 2 + 1


def sigma_rob(data, iterations=1, thresh=3.0, axis=None):
    """
    Iterative m.a.d. based sigma with positive outlier rejection.
    """
    noise = mad_std(data, axis=axis)
    for _ in range(iterations):
        ind = (np.abs(data) <= thresh * noise).nonzero()
        noise = mad_std(data[ind], axis=axis)
    return noise


cube = SpectralCube.read(fourteenB_HI_data_path(cube_name))

# noise = Noise(cube)

# scale = noise.scale

# scale = sig_clip(cube[-1].value, nsig=10)

pixscale = proj_plane_pixel_scales(cube.wcs)[0]

# cube = cube.with_mask(cube > 3 * scale * u.Jy)

# # Want to smooth the mask edges
mask = cube.mask.include()

# Perform smooth and mask on spectra
posns = np.where(mask.sum(0) > 0)

# Set smoothing parameters and # consecutive channels.
# Real emission has fairly broad profiles here. Smooth over 6 km/s.
smooth_chans = int(round_up_to_odd(6 / 0.2))
# 10 consecutive channels must be above the MAD level to be real emission.
num_chans = 30
# Smoothed spectrum must have a peak intensity of >=3 * MAD to be considered
# real.
peak_snr = 5
# Cutoff level
min_snr = 2
# Where to cut at the line edges
edge_thresh = 1

# In case single spectra need to be inspected.
verbose = False

for i, j in ProgressBar(zip(*posns)):

    spectrum = cube[:, i, j].value

    # Assume for now that the noise level is ~ constant across the channels.
    # Not quite true, but good enough just for creating a mask.
    smoothed = medfilt(spectrum, smooth_chans)

    mad = sigma_rob(smoothed, thresh=min_snr, iterations=5)

    good_posns = np.where(smoothed > min_snr * mad)[0]

    if good_posns.size < num_chans:
        mask[:, i, j] = False
        continue

    sequences = []
    for k, g in groupby(enumerate(good_posns), lambda (i, x): i - x):
        sequences.append(map(itemgetter(1), g))

    # Check length and peak
    sequences = [seq for seq in sequences if len(seq) >= num_chans and
                 np.nanmax(smoothed[seq]) >= peak_snr * mad]

    if len(sequences) == 0:
        mask[:, i, j] = False
        continue

    # Now take each valid sequence and expand the edges until the smoothed
    # spectrum approaches zero.
    edges = [[seq[0], seq[-1]] for seq in sequences]
    for n, edge in enumerate(edges):
        # Lower side
        if n == 0:
            start_posn = edge[0]
            stop_posn = 0
        else:
            start_posn = edge[0] - edges[n - 1][0]
            stop_posn = edges[n - 1][0]

        for pt in np.arange(start_posn, stop_posn, -1):
            if smoothed[pt] <= mad * edge_thresh:
                break

            sequences[n].insert(0, pt)

        # Upper side
        start_posn = edge[1]
        if n == len(edges) - 1:
            stop_posn = cube.shape[0]
        else:
            stop_posn = edges[n + 1][0]

        for pt in np.arange(start_posn, stop_posn, 1):
            if smoothed[pt] <= mad * edge_thresh:
                break

            sequences[n].insert(0, pt)

    bad_posns = \
        list(set(np.arange(cube.shape[0])) - set(list(chain(*sequences))))

    mask[:, i, j][bad_posns] = False

    if verbose:
        p.plot(spectrum)
        p.plot(smoothed)
        p.vlines(np.where(mask[:, i, j])[0][-1], 0, np.nanmax(spectrum))
        p.vlines(np.where(mask[:, i, j])[0][0], 0, np.nanmax(spectrum))
        p.plot(smoothed * mask[:, i, j], 'bD')
        p.draw()
        raw_input("Next spectrum?")

kernel = average_beams(cube.beams).as_tophat_kernel(pixscale)
kernel_pix = (kernel.array > 0).sum()

for i in ProgressBar(mask.shape[0]):
    mask[i] = nd.binary_opening(mask[i], kernel)
    mask[i] = nd.binary_closing(mask[i], kernel)
    mask[i] = mo.remove_small_objects(mask[i], min_size=kernel_pix,
                                      connectivity=2)


new_header = cube.header.copy()
new_header["BUNIT"] = ""
new_header["BITPIX"] = 8

mask_hdu = fits.PrimaryHDU(mask.astype('>i2'), header=new_header)
mask_hdu.writeto(fourteenB_HI_data_path(mask_name, no_check=True),
                 clobber=True)

# print("Now the source mask for the rotation subtracted cube.")
# cube = SpectralCube.read(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits"))
# scale = sig_clip(cube[-1].value, nsig=10)

# pixscale = proj_plane_pixel_scales(cube.wcs)[0]

# cube = cube.with_mask(cube > 5 * scale * u.Jy)

# # Want to smooth the mask edges
# mask = cube.mask.include()

# kernel = average_beams(cube.beams).as_tophat_kernel(pixscale)

# for i in ProgressBar(xrange(mask.shape[0])):
#     mask[i] = nd.binary_opening(mask[i], kernel)
#     mask[i] = nd.binary_closing(mask[i], kernel)

#     mask[i] = nd.binary_dilation(mask[i], mo.disk(10))

#     # Add peak brightness in region check to remove spurious small features

# new_header = cube.header.copy()
# new_header["BUNIT"] = ""
# new_header["BITPIX"] = 8

# mask_hdu = fits.PrimaryHDU(mask.astype('>i2'), header=new_header)
# mask_hdu.writeto(fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub_source_mask.fits",
#                                         no_check=True),
#                  clobber=True)
