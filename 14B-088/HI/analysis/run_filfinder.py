
from fil_finder import fil_finder_2D
from fil_finder.width import gauss_model
from basics import BubbleFinder2D
from spectral_cube.lower_dimensional_structures import Projection
from astropy.io import fits
from radio_beam import Beam
from astropy.wcs import WCS
import astropy.units as u
import matplotlib.pyplot as p
import scipy.ndimage as nd
from scipy.stats import binned_statistic
import numpy as np
from skimage.morphology import medial_axis
from galaxies import Galaxy
from warnings import warn


from HI_veldisp_profile import radial_profile

'''
Filaments in M33? Why not?
'''


def fit_skeleton_width(skel_array, intensity_array, pixscale, beam_width,
                       max_dist=50):
    '''
    Reproduce the width fitting of FilFinder, but simplified for a single
    super-structure.

    Parameters
    ----------
    pixscale : float
        Conversion of pixels to pc.
    beam_width : float
        Beam size in pc.
    max_dist : float or int
        Maximum radial distance in PIXELS.
    '''

    dist_arr = nd.distance_transform_edt(~skel_array)

    # bins = np.linspace(dist_arr.min(), max_dist,
    #                    int(np.sqrt(np.isfinite(dist_arr).sum())))

    del_bin = 0.25 * (beam_width / pixscale)

    bins = np.arange(dist_arr.min(), max_dist + del_bin,
                     del_bin)

    vals, bin_edges, bin_num = \
        binned_statistic(dist_arr.ravel(), intensity_array.ravel(),
                         bins=bins,
                         statistic=np.nanmean)

    variances = \
        binned_statistic(dist_arr.ravel(), intensity_array.ravel(),
                         bins=bins,
                         statistic=np.nanvar)[0]

    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width / 2

    # weights = np.array([(bin_num == i).sum() for i in range(len(bins))])
    # weights = weights[weights > 0]

    # Convert the bin_centers to pc
    bin_centers *= pixscale

    # Nans are returned when nothing usable in a bin. Remove these before
    # fitting
    bin_centers = bin_centers[np.isfinite(vals)]
    variances = variances[np.isfinite(vals)]
    vals = vals[np.isfinite(vals)]

    weights = 1 / variances

    profile = np.array([bin_centers, vals, np.sqrt(variances)])

    # Now fit to a gaussian
    fit, err, model, names, fail = \
        gauss_model(bin_centers, vals, weights, beam_width)

    if fail:
        warn("Fail flag was raised. Check the output")
        # raise ValueError("Fit failed.")

    return fit, err, model, profile


mom0_fits = fits.open("/home/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits")[0]
beam = Beam.from_fits_header(mom0_fits.header)
mom0 = Projection(mom0_fits.data * beam.jtok(1.414 * u.GHz) / 1000. * u.km / u.s,
                  wcs=WCS(mom0_fits.header))
mom0.meta['beam'] = beam

# Create the bubble mask instead of letting FilFinder to do it.
bub = BubbleFinder2D(mom0, sigma=80. * beam.jtok(1.414 * u.GHz) / 1000.)

# fils = fil_finder_2D(mom0.value, mom0.header, 10, distance=0.84e6)
# fils.mask = ~(bub.mask.copy())
# fils.medskel()
# fils.analyze_skeletons()
# # So at least on of the radial profiles fails. BUT the second fit is to a
# # skeleton that is essentially the entire disk, so plot without interactivity
# # and save the plot and the parameters shown in verbose mode.
# p.ioff()
# fils.find_widths(verbose=True, max_distance=500, auto_cut=False, try_nonparam=False)

# Fit Parameters: [ 541.31726502  129.85351117  180.0710914   304.01262168
# Fit Errors: [ 0.89151974  0.48394493  0.27313627  1.1462345 ]

skeleton = medial_axis(~bub.mask)

# Overplot the skeleton on the moment0
ax = p.subplot(111, projection=mom0.wcs)
ax.imshow(mom0.value, origin='lower')
ax.contour(skeleton, colors='r')
p.draw()

raw_input("Next plot?")
p.clf()

pixscale = \
    mom0.header['CDELT2'] * (np.pi / 180.) * 0.84e6  # * u.pc
beam_width = \
    mom0.meta['beam'].major.to(u.deg).value * (np.pi / 180.) * 0.84e6  # * u.pc

# HI profile and fit
fit, err, model, profile = \
    fit_skeleton_width(skeleton, mom0.value, pixscale, beam_width)

rads = np.linspace(0, profile[0].max(), 100)

p.plot(profile[0], profile[1], 'bD')
# p.errorbar(profile[0], profile[1], yerr=profile[2], fmt="D",
#            color="b")
p.plot(rads, model(rads, *fit), 'r-')
p.xlabel("Radial Distance (pc)")
p.ylabel("Integrated Intensity (K km/s)")
p.grid()
p.draw()

raw_input("Next plot?")
p.clf()


# Now use the regridded integrated CO map with the same skeleton
# Created with co_comparison/co_reproject.py
mom0_co_fits = fits.open("/media/eric/Data_3/M33/co21/m33.ico.hireprojection.fits")[0]
# Correction for beam efficiency added in.
mom0_co = Projection(mom0_co_fits.data / 0.75, wcs=WCS(mom0_co_fits.header))
mom0_co.meta['beam'] = Beam.from_fits_header(mom0_co_fits.header)

# On the same grid, so pixscales are the same
beam_width_co = \
    mom0_co.meta['beam'].major.to(u.deg).value * (np.pi / 180.) * 0.84e6  # * u.pc

fit_co, err_co, model_co, profile_co = \
    fit_skeleton_width(skeleton, mom0_co.value, pixscale, beam_width_co)

# p.errorbar(profile_co[0], profile_co[1], yerr=profile_co[2], fmt="D",
#            color="b")
p.plot(profile_co[0], profile_co[1], "bD")
p.plot(rads, model(rads, *fit_co), 'r-')
p.xlabel("Radial Distance (pc)")
p.ylabel("Integrated Intensity (K km/s)")
p.grid()
p.draw()

raw_input("Next plot?")
p.clf()

# fils = fil_finder_2D(mom0_co.value, mom0.header, 10, distance=0.84e6)
# fils.mask = ~(bub.mask.copy())
# fils.medskel()
# fils.analyze_skeletons()
# fils.find_widths(verbose=True, max_distance=500, auto_cut=False, try_nonparam=False)

# Now scale the HI and CO to have the same amplitude and ignore the HI bkg

p.plot(profile[0], model(profile[0], *[1, fit[1], 0]), 'b-', label="HI Model")
p.plot(profile[0], model(profile[0], *[1, fit_co[1], 0]), 'g--',
       label='CO Model')
p.legend()
p.grid()
p.draw()

raw_input("Next plot?")
p.clf()

# Now create a radial profile of the distance between filaments
g = Galaxy("M33")

# Create a radial profile of the HI vel disp out to 8 kpc.
dist_arr = Projection(nd.distance_transform_edt(~skeleton) * pixscale,
                      unit=u.pc, wcs=mom0.wcs)
dist_arr.meta['beam'] = beam

rs, sd, sd_sigma = radial_profile(g, dist_arr, max_rad=10 * u.kpc,
                                  dr=100 * u.pc)

# Multiply 2 to get the average distance instead of radius
sd = 2 * sd
sd_sigma = 2 * sd_sigma

p.errorbar(rs.value, sd.value,
           yerr=sd_sigma.value, fmt="-", color="b",
           drawstyle='steps-mid')
p.xlabel("R (kpc)")
p.ylabel("Average Distance between Filaments (pc)")
p.grid()
p.xlim([0, 8])
p.ylim([80, 300])
p.draw()
