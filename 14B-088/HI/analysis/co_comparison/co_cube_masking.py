
'''
Apply the signal masking procedure on the IRAM data, using the spatial RMS
map
'''

import astropy.units as u
from astropy.convolution import Box1DKernel
from signal_id import Noise
from scipy import ndimage as nd
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.utils.console import ProgressBar
import skimage.morphology as mo
import numpy as np
from radio_beam import Beam
from itertools import groupby, chain
from operator import itemgetter
import matplotlib.pyplot as p


def round_up_to_odd(f):
    return np.ceil(f) // 2 * 2 + 1


def make_signal_mask(cube, smooth_chans=3, min_chan=7, peak_snr=5.,
                     min_snr=3.5, edge_thresh=1.5, verbose=False,
                     noise_map=None):
    '''
    Create a robust signal mask by requiring spatial and spectral
    connectivity.
    '''

    pixscale = proj_plane_pixel_scales(cube.wcs)[0]

    # # Want to smooth the mask edges
    mask = cube.mask.include().copy()

    # Set smoothing parameters and # consecutive channels.
    smooth_chans = int(round_up_to_odd(smooth_chans))

    # consecutive channels to be real emission.
    num_chans = min_chan

    # Smooth the cube, then create a noise model
    spec_kernel = Box1DKernel(smooth_chans)
    smooth_cube = cube.spectral_smooth(spec_kernel)

    if noise_map is not None:
        noise = Noise(smooth_cube)
        noise.estimate_noise(spectral_flat=True)
        noise.get_scale_cube()

        snr = noise.snr.copy()

    else:
        try:
            snr = cube.filled_data[:] / noise_map
        except Exception as e:
            print(e)
            print("You may have to allow for huge cube operations.")

    snr[np.isnan(snr)] = 0.0

    posns = np.where(snr.max(axis=0) >= min_snr)

    # Blank the spectra for which none are above min_snr
    bad_pos = np.where(snr.max(axis=0) < min_snr)
    mask[:, bad_pos[0], bad_pos[1]] = False

    for i, j in ProgressBar(zip(*posns)):

        # Look for all pixels above min_snr
        good_posns = np.where(snr[:, i, j] > min_snr)[0]

        # Reject if the total is less than connectivity requirement
        if good_posns.size < num_chans:
            mask[:, i, j] = False
            continue

        # Find connected pixels
        sequences = []
        for k, g in groupby(enumerate(good_posns), lambda (i, x): i - x):
            sequences.append(map(itemgetter(1), g))

        # Check length and peak. Require a minimum of 3 pixels above the noise
        # to grow from.
        sequences = [seq for seq in sequences if len(seq) >= 3 and
                     np.nanmax(snr[:, i, j][seq]) >= peak_snr]

        # Continue if no good sequences found
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
                # if smoothed[pt] <= mad * edge_thresh:
                if snr[:, i, j][pt] <= edge_thresh:
                    break

                sequences[n].insert(0, pt)

            # Upper side
            start_posn = edge[1]
            if n == len(edges) - 1:
                stop_posn = cube.shape[0]
            else:
                stop_posn = edges[n + 1][0]

            for pt in np.arange(start_posn, stop_posn, 1):
                # if smoothed[pt] <= mad * edge_thresh:
                if snr[:, i, j][pt] <= edge_thresh:
                    break

                sequences[n].insert(0, pt)

        # Final check for the min peak level and ensure all meet the
        # spectral connectivity requirement
        sequences = [seq for seq in sequences if len(seq) >= num_chans and
                     np.nanmax(snr[:, i, j][seq]) >= peak_snr]

        if len(sequences) == 0:
            mask[:, i, j] = False
            continue

        bad_posns = \
            list(set(np.arange(cube.shape[0])) - set(list(chain(*sequences))))

        mask[:, i, j][bad_posns] = False

        if verbose:
            p.subplot(121)
            p.plot(cube.spectral_axis.value, snr[:, i, j])
            min_val = cube.spectral_axis.value[np.where(mask[:, i, j])[0][-1]]
            max_val = cube.spectral_axis.value[np.where(mask[:, i, j])[0][0]]
            p.vlines(min_val, 0,
                     np.nanmax(snr[:, i, j]))
            p.vlines(max_val, 0,
                     np.nanmax(snr[:, i, j]))
            p.plot(cube.spectral_axis.value,
                   snr[:, i, j] * mask[:, i, j], 'bD')

            p.subplot(122)
            p.plot(cube.spectral_axis.value, cube[:, i, j], label='Cube')
            p.plot(cube.spectral_axis.value, smooth_cube[:, i, j],
                   label='Smooth Cube')
            p.axvline(min_val)
            p.axvline(max_val)
            p.plot(cube.spectral_axis.value,
                   smooth_cube[:, i, j] * mask[:, i, j], 'bD')
            p.draw()
            raw_input("Next spectrum?")
            p.clf()

    # initial_mask = mask.copy()

    # Now set the spatial connectivity requirements.

    kernel = cube.beam.as_tophat_kernel(pixscale)
    # kernel = Beam(major=0.75 * cube.beam.major, minor=0.75 * cube.beam.minor,
    #               pa=cube.beam.pa).as_tophat_kernel(pixscale)
    kernel_pix = (kernel.array > 0).sum()

    # Avoid edge effects in closing by padding by 1 in each axis
    mask = np.pad(mask, ((0, 0), (1, 1), (1, 1)), 'constant',
                  constant_values=False)

    for i in ProgressBar(mask.shape[0]):
        mask[i] = nd.binary_opening(mask[i], kernel)
        mask[i] = nd.binary_closing(mask[i], kernel)
        mask[i] = mo.remove_small_objects(mask[i], min_size=kernel_pix,
                                          connectivity=2)
        mask[i] = mo.remove_small_holes(mask[i], min_size=kernel_pix,
                                        connectivity=2)

    # Remove padding
    mask = mask[:, 1:-1, 1:-1]

    # Each region must contain a point above the peak_snr
    labels, num = nd.label(mask, np.ones((3, 3, 3)))
    for n in range(1, num + 1):
        pts = np.where(labels == n)
        if np.nanmax(snr[pts]) < peak_snr:
            mask[pts] = False

    masked_cube = cube.with_mask(mask)

    return masked_cube


def simple_masking(cube, rms, min_sig=3, max_sig=5, min_pix=27):
    '''
    Find connected regions above 3 sigma that contain a pixel at least above
    5 sigma, and contains some minimum number of pixels.
    '''

    mask_low = (cube > min_sig * rms * cube.unit).include()
    mask_high = (cube > max_sig * rms * cube.unit).include()

    mask_low = mo.remove_small_objects(mask_low, min_size=min_pix,
                                       connectivity=2)

    # Remove all regions that do not contain a 5 sigma pixel
    kernel = np.ones((3, 3, 3))
    labels, num = nd.label(mask_low, kernel)

    for i in xrange(1, num + 1):
        pix = np.where(labels == i)
        if np.any(mask_high[pix]):
            continue
        mask_low[pix] = False

    return cube.with_mask(mask_low)


if __name__ == "__main__":

    from spectral_cube import SpectralCube
    import astropy.io.fits as fits
    from reproject import reproject_interp

    from paths import (iram_co21_data_path, fourteenB_HI_data_path,
                       paper1_figures_path)
    from constants import moment0_name


    cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
    rms = fits.open(iram_co21_data_path("m33.rms.fits"))[0].data.squeeze()
    rms[rms < 0] = np.NaN

    masked_cube = simple_masking(cube, rms)
    masked_cube.write(iram_co21_data_path("m33.co21_iram.masked.fits", no_check=True))

    mom0 = masked_cube.moment0()
    mom0.write(iram_co21_data_path("m33.co21_iram.masked.moment0.fits", no_check=True))

    mom1 = masked_cube.moment1()
    mom1.write(iram_co21_data_path("m33.co21_iram.masked.moment1.fits", no_check=True))

    lwidth = masked_cube.linewidth_sigma()
    lwidth.write(iram_co21_data_path("m33.co21_iram.masked.lwidth.fits", no_check=True))

    # Make HI reprojected versions
    hi_mom0 = fits.open(fourteenB_HI_data_path(moment0_name))[0]
    beam = Beam.from_fits_header(hi_mom0.header)

    masked_cube_conv = masked_cube.convolve_to(beam)
    mom0 = masked_cube_conv.moment0()
    mom1 = masked_cube_conv.moment1()
    lwidth = masked_cube_conv.linewidth_sigma()

    mom0_reproj = reproject_interp(mom0.hdu, hi_mom0.header)[0]
    mom0_reproj_hdu = fits.PrimaryHDU(mom0_reproj, hi_mom0.header)
    mom0_reproj_hdu.writeto(iram_co21_data_path("m33.co21_iram.masked.moment0.hireprojection.fits", no_check=True))

    mom1_reproj = reproject_interp(mom1.hdu, hi_mom0.header)[0]
    mom1_reproj_hdu = fits.PrimaryHDU(mom1_reproj, hi_mom0.header)
    mom1_reproj_hdu.writeto(iram_co21_data_path("m33.co21_iram.masked.moment1.hireprojection.fits", no_check=True))

    lwidth_reproj = reproject_interp(lwidth.hdu, hi_mom0.header)[0]
    lwidth_reproj_hdu = fits.PrimaryHDU(lwidth_reproj, hi_mom0.header)
    lwidth_reproj_hdu.writeto(iram_co21_data_path("m33.co21_iram.masked.lwidth.hireprojection.fits", no_check=True))
