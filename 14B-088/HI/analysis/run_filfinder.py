
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
from warnings import warn


from HI_veldisp_profile import radial_profile
from constants import hi_freq

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


if __name__ == "__main__":

    from analysis.paths import (fourteenB_HI_data_path, paper1_figures_path,
                                iram_co21_data_path)
    from analysis.constants import moment0_name
    from analysis.galaxy_params import gal

    mom0_fits = fits.open(fourteenB_HI_data_path(moment0_name))[0]
    beam = Beam.from_fits_header(mom0_fits.header)
    mom0 = Projection(mom0_fits.data * beam.jtok(hi_freq) / 1000. *
                      u.km / u.s,
                      wcs=WCS(mom0_fits.header))
    mom0.meta['beam'] = beam

    # Create the bubble mask instead of letting FilFinder to do it.
    bub = BubbleFinder2D(mom0, sigma=60)
    bub.create_mask(bkg_nsig=6, region_min_nsig=7, fill_radius=1)

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

    # Reject regions that don't have enough pixels
    beam_pix = int(np.ceil(beam.major.value / mom0.header["CDELT2"]))
    labels, n = nd.label(skeleton, np.ones((3, 3)))
    num_pixels = nd.sum(skeleton, labels, range(1, n + 1))
    for i, pix in enumerate(num_pixels):
        # Longer than ~3 beam widths.
        if pix < 3 * beam_pix:
            skeleton[np.where(labels == i + 1)] = False

    # Overplot the skeleton on the moment0
    ax = p.subplot(111, projection=mom0.wcs)
    ax.imshow(mom0.value, origin='lower')
    ax.contour(skeleton, colors='r')
    p.draw()

    p.savefig(paper1_figures_path("moment0_w_skeletons.pdf"), rasterize=True)
    p.savefig(paper1_figures_path("moment0_w_skeletons.png"))

    # raw_input("Next plot?")
    p.clf()

    pixscale = \
        mom0.header['CDELT2'] * (np.pi / 180.) * gal.distance.to(u.pc).value
    beam_width = \
        mom0.meta['beam'].major.to(u.deg).value * (np.pi / 180.) * \
        gal.distance.to(u.pc).value

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

    p.savefig(paper1_figures_path("filfinder_radial_profile_hi.pdf"))
    p.savefig(paper1_figures_path("filfinder_radial_profile_hi.png"))

    fraction_in_peak = np.sqrt(2 * np.pi) * (fit[0] - fit[2]) * \
        (fit[1] / pixscale) * skeleton.sum() / np.nansum(mom0.value)

    print("Fraction of emission above the background level: {}"
          .format(fraction_in_peak))
    # 0.543

    # raw_input("Next plot?")
    p.clf()

    # Now use the regridded integrated CO map with the same skeleton
    # Created with co_comparison/co_reproject.py
    mom0_co_fits = \
        fits.open(iram_co21_data_path("m33.ico.hireprojection.fits"))[0]
    # Correction for beam efficiency added in.
    mom0_co = \
        Projection(mom0_co_fits.data / 0.75, wcs=WCS(mom0_co_fits.header))
    mom0_co.meta['beam'] = Beam.from_fits_header(mom0_co_fits.header)

    # On the same grid, so pixscales are the same
    beam_width_co = \
        mom0_co.meta['beam'].major.to(u.deg).value * (np.pi / 180.) * \
        gal.distance.to(u.pc).value

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

    p.savefig(paper1_figures_path("filfinder_radial_profile_co.pdf"))
    p.savefig(paper1_figures_path("filfinder_radial_profile_co.png"))

    fraction_in_peak_co = np.sqrt(2 * np.pi) * (fit_co[0]) * \
        (fit_co[1] / pixscale) * skeleton.sum() / np.nansum(mom0_co.value)

    print("Fraction of CO emission above the background level: {}"
          .format(fraction_in_peak_co))
    # 1.12 Yeah, this one doesn't work. And it doesn't matter since this should
    # pretty much be all of it.

    # raw_input("Next plot?")
    p.clf()

    # fils = fil_finder_2D(mom0_co.value, mom0.header, 10, distance=0.84e6)
    # fils.mask = ~(bub.mask.copy())
    # fils.medskel()
    # fils.analyze_skeletons()
    # fils.find_widths(verbose=True, max_distance=500, auto_cut=False, try_nonparam=False)

    # Now scale the HI and CO to have the same amplitude and ignore the HI bkg

    p.plot(profile[0], model(profile[0], *[1, fit[1], 0]), 'b-',
           label="HI Model")
    p.plot(profile[0], model(profile[0], *[1, fit_co[1], 0]), 'g--',
           label='CO Model')
    p.xlabel("Radial Distance (pc)")
    p.ylabel("Normalized Intensity")
    p.legend()
    p.grid()
    p.draw()

    # raw_input("Next plot?")
    p.savefig(paper1_figures_path("filfinder_radial_profile_models.pdf"))
    p.savefig(paper1_figures_path("filfinder_radial_profile_models.png"))

    p.clf()

    # Now create a radial profile of the distance between filaments

    # Create a radial profile of the HI vel disp out to 8 kpc.
    dist_arr = Projection(nd.distance_transform_edt(~skeleton) * pixscale,
                          unit=u.pc, wcs=mom0.wcs)
    dist_arr.meta['beam'] = beam

    rs, sd, sd_sigma = radial_profile(gal, dist_arr, max_rad=10 * u.kpc,
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

    p.savefig(paper1_figures_path("filfinder_skeleton_distance_radial_profile.pdf"))
    p.savefig(paper1_figures_path("filfinder_skeleton_distance_radial_profile.png"))

    p.close()
