
'''
Highlight the peaks from the NOEMA arm map on the HI.

Throw on some A config continuum peaks for good measure.
'''

from spectral_cube import SpectralCube
import astropy.io.fits as fits
import numpy as np
import scipy.ndimage as nd
from basics.utils import sig_clip
import astropy.units as u
from aplpy import FITSFigure
import os


data_path = "/mnt/MyRAID/M33/VLA/14B-088/combined_HI/"
cube = \
    SpectralCube.read(os.path.join(data_path,
                                   "M33_14B-088_AT0206_HI.clean.image.fits"))

# Lots of crap on the edge since the B-config fields don't cover all of the C
pbcov = fits.open(os.path.join(data_path,
                               "M33_HI_14B-088_AT0206_pbcov.fits"))[0].data[0]
cube = cube.with_mask(pbcov > 0.6).minimal_subcube()

# This channel is nearly empty of signal.
sigma = sig_clip(cube[-1].value, nsig=10)

# Lazy zeroth moment
mom0 = cube.with_mask(cube > 2 * sigma * u.Jy).moment0()

# fig = FITSFigure("/mnt/MyRAID/M33/VLA/14B-088/HI/full_imaging_noSD/M33_14B-088_HI.clean.image.pbcov_gt_0.5_masked.mom0.fits")
fig = FITSFigure(mom0.hdu)

fig.show_grayscale()

# The NOEMA map has not (yet) been cleaned, so just clip to the tops of the
# brightest peaks! Specifically, highlight that super-duper bright clump
co_cube = SpectralCube.read("/home/eric/bigdata/ekoch/M33/co21_noema/M33-ARM13-merged.fits")
co21_mom0 = co_cube.with_mask(co_cube > 0.3 * u.Jy).moment0()

# And an unmasked so we can show the map coverage.
co21_mom0_mask = co_cube.moment0().value > 0.005
co21_mom0_mask = nd.binary_fill_holes(co21_mom0_mask)
co21_mom0_mask = nd.binary_closing(co21_mom0_mask, np.ones((12, 12)))


fig.show_contour(co21_mom0.hdu)
fig.show_contour(fits.PrimaryHDU(co21_mom0_mask.astype(int), header=co21_mom0.header))

fig.save(os.path.expanduser("~/Dropbox/Various Plots/Proposals/noema_dirty_on_at0206_14b088_hi.pdf"))
fig.save(os.path.expanduser("~/Dropbox/Various Plots/Proposals/noema_dirty_on_at0206_14b088_hi.png"))
