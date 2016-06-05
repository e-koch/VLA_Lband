
from fil_finder import fil_finder_2D
from basics import BubbleFinder2D
from spectral_cube.lower_dimensional_structures import Projection
from astropy.io import fits
from radio_beam import Beam
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as p

'''
Filaments in M33? Why not?
'''

mom0_fits = fits.open("/home/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits")[0]
mom0 = Projection(mom0_fits.data, wcs=WCS(mom0_fits.header))
mom0.meta['beam'] = Beam.from_fits_header(mom0_fits.header)

# Create the bubble mask instead of letting FilFinder to do it.
bub = BubbleFinder2D(mom0, sigma=80.)

fils = fil_finder_2D(mom0.value, mom0.header, 10, distance=0.84e6)
fils.mask = ~(bub.mask.copy())
fils.medskel()
fils.analyze_skeletons()
# So at least on of the radial profiles fails. BUT the second fit is to a
# skeleton that is essentially the entire disk, so plot without interactivity
# and save the plot and the parameters shown in verbose mode.
p.ioff()
fils.find_widths(verbose=True, max_distance=500, auto_cut=False, try_nonparam=False)

# Fit Parameters: [ 541.31726502  129.85351117  180.0710914   304.01262168
# Fit Errors: [ 0.89151974  0.48394493  0.27313627  1.1462345 ]
