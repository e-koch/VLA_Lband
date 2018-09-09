
'''
Pull out HI properties (and/or others) from a set of point sources.

Create a distance map as a function of distance from the nearest source.
'''

import astropy.coordinates as coord
from astropy.table import Table, Column
import astropy.units as u
import astropy.constants as const
import numpy as np
from galaxies import Galaxy
import scipy.ndimage as nd
from astropy.io import fits
from spectral_cube import SpectralCube, Projection
from spectral_cube.analysis_utilities import stack_spectra
from astropy.utils.console import ProgressBar
import matplotlib.pyplot as plt

from plotting_styles import default_figure
from constants import hi_freq, hi_mass_conversion
# from paths import allfigs_path


def distance_map_from_catalogue(gal, tab, header, ra_key="RA", dec_key="Dec",
                                diam_key=None):
    '''
    Create a distance map from a set of sky location in a catalogue.
    '''

    if not isinstance(gal, Galaxy):
        raise TypeError("gal must be a Galaxy instance.")

    ra = tab[ra_key]
    dec = tab[dec_key]

    coords = coord.SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg))

    # Assumes that the table column has a unit attached that Table can distinguish
    if diam_key is not None:
        # Assume pc. Lost units in the save table??
        diams = tab[diam_key].quantity * u.pc

        radii = gal.radius(header=header)

    coord_map = gal.skycoord_grid(header=header)

    object_mask = np.zeros_like(coord_map.ra.value, dtype=int)

    # Loop through and mask points belonging at a remnant, or the nearest point
    for i, co in enumerate(coords):

        mask_index = np.unravel_index(coord_map.separation(co).argmin(),
                                      object_mask.shape)

        if diam_key is not None:

            # Major axis diameter

            diam_rad = (diams[i].to(u.pc) / gal.distance).to(u.dimensionless_unscaled).value * u.rad
            diam_pix = diam_rad.to(u.deg).value / np.abs(header['CDELT2'])

            # Gather all pixels with a circular region
            yy, xx = np.mgrid[-(int(diam_pix)//2 + 1):int(diam_pix)//2 + 1,
                              -(int(diam_pix)//2 + 1):int(diam_pix)//2 + 1]

            # Find all pixels within the diameter
            valids = np.where(np.sqrt(yy**2 + xx**2) < diam_pix / 2.)

            y_pts = valids[0] + mask_index[0]
            x_pts = valids[1] + mask_index[1]

            mask_index = (y_pts, x_pts)

        object_mask[mask_index] = i + 1
        # print(object_mask[mask_index])
        # print(mask_index)
        # print((object_mask > 0).sum())

    dist_transf = nd.distance_transform_edt(~(object_mask > 0))

    return object_mask, dist_transf


def find_bubble_props(dist_bins, int_profile, lwidth_profile, obj_diam,
                      disk_height=100 * u.pc / np.cos(55.1 * u.deg),
                      mass_conv_factor=None):
    '''
    Dumb estimations of bubble properties based on integrated intensity and
    line width profiles.
    '''

    # Define the shell radius based on the distance of the peak

    arg_max = np.argmax(int_profile)

    # If the centre is the peak, assume it is unresolved
    if arg_max == 0:
        shell_rad = obj_diam / 2.
    else:
        shell_rad = obj_diam / 2. + dist_bins[arg_max]

    # Assume a disk scale height and check if the radius of the shell
    # exceeds it
    if shell_rad > disk_height:
        # It has maybe broken out of the disk. Adjust volume as needed
        # Subtract off caps of the sphere
        vol = (4 * np.pi / 3.) * shell_rad**3 - \
            (2 * np.pi / 3.) * (shell_rad - disk_height)**2 * (2 * shell_rad + disk_height)

    else:
        # Likely still contained within the disk
        vol = (4 * np.pi / 3.) * shell_rad**3

    # Awful estimations of the velocity expansion. Assume velocity dispersion
    # is exactly the same...
    # Don't know how to do that with any sort of logic applied, so let it be
    # the dispersion in the peak bin.
    v_exp = lwidth_profile[arg_max]

    # Now the integ intensity. If unresolved, we don't have an estimate of the
    # background. Assume the last distance bin as a background?? Otherwise take
    # the larger of the innermost and outermost when resolved.

    peak_int = int_profile[arg_max]

    if arg_max == 0:
        bkg_int = int_profile[-1]
    else:
        bkg_int = max(int_profile[0], int_profile[-1])

    hole_mass = np.pi * shell_rad**2 * bkg_int

    shell_mass = np.pi * shell_rad**2 * \
        (peak_int - bkg_int)

    if mass_conv_factor is not None:
        hole_mass *= mass_conv_factor
        shell_mass *= mass_conv_factor

    # Estimate an avg volume density within the hole. Don't do this
    # for unresolved holes
    if arg_max == 0:
        energy = np.NaN * u.erg
        vol_dens = np.NaN * u.cm**-3
    else:
        # Chevalier 74 expansion energy formula
        vol_dens = ((shell_mass / (1.4 * const.m_p)) / vol).to(u.cm**-3)
        energy = 5.3e43 * vol_dens.value**1.12 * \
            shell_rad.to(u.pc).value**3.12 * v_exp.to(u.km / u.s).value**1.4 * u.erg

    return shell_rad, vol, v_exp, hole_mass, shell_mass, vol_dens, energy

default_figure()

# Update this for the server files (geometry should be the same though)
gal = Galaxy("M33")
gal.distance = 840 * u.kpc

hi_cube = SpectralCube.read("/Volumes/Travel_Data/M33_2/HI/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.fits")

peak_vel = Projection.from_hdu(fits.open("/Volumes/Travel_Data/M33_2/HI/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peakvels.fits"))
mom0 = Projection.from_hdu(fits.open("/Volumes/Travel_Data/M33_2/HI/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.mom0.fits"))

beam = mom0.beam
moment0_Kkm_s = beam.jtok(hi_freq).value * mom0.value / 1000.
moment0_surfdens = moment0_Kkm_s * hi_mass_conversion * (u.K * u.km / u.s) * np.cos(55.1 * u.deg)


lwidth = Projection.from_hdu(fits.open("/Volumes/Travel_Data/M33_2/HI/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.lwidth.fits"))

snr_tab = Table.read("/Volumes/Travel_Data/M33_2/MMT_SNR_catalogue_long18_combined.txt",
                     format='ascii')

# Also consider weighting by something like ~1/sqrt(L) to place distances
# on a common "scale"
index_mask, dist_transf = \
    distance_map_from_catalogue(gal, snr_tab, hi_cube.header,
                                diam_key='D')

# Get all points within ~100 pc.

dist_limit = np.arange(10) * 100 * u.pc

stacked_spectra = []
lwidth_bins = []
intint_bins = []

# Pick out individual regions
num = index_mask.max()

for n in ProgressBar(range(1, num + 1)):

    reg_mask = index_mask == n

    dist_transf_reg = nd.distance_transform_edt(~reg_mask)

    lwidth_reg = []
    intint_reg = []

    # Calculate avg properties within the region
    lwidth_reg.append([np.nanmean(lwidth[reg_mask].value),
                       np.nanstd(lwidth[reg_mask].value) / np.sqrt(reg_mask.sum() / 41.)])
    intint_reg.append([np.nanmean(moment0_surfdens[reg_mask].value),
                       np.nanstd(moment0_surfdens[reg_mask].value) / np.sqrt(reg_mask.sum() / 41.)])

    for i, (low, high) in enumerate(zip(dist_limit[:-1], dist_limit[1:])):

        # print("On bin {}".format(i + 1))

        dist_ang_low = (low / gal.distance.to(u.pc)).value * u.rad
        dist_pix_low = dist_ang_low.to(u.deg).value / np.abs(hi_cube.header['CDELT2'])

        dist_ang_high = (high / gal.distance.to(u.pc)).value * u.rad
        dist_pix_high = dist_ang_high.to(u.deg).value / np.abs(hi_cube.header['CDELT2'])

        dist_mask = np.logical_and(dist_transf_reg > dist_pix_low,
                                   dist_transf_reg <= dist_pix_high)

        num_beams = dist_mask.sum() / 41.

        intint_reg.append([np.nanmean(moment0_surfdens[dist_mask].value),
                           np.nanstd(moment0_surfdens[dist_mask].value) / np.sqrt(num_beams)])
        lwidth_reg.append([np.nanmean(lwidth[dist_mask].value),
                           np.nanstd(lwidth[dist_mask].value) / np.sqrt(num_beams)])

        # stacked_spectra.append(stack_spectra(hi_cube, peak_vel,
        #                                      xy_posns=np.where(dist_mask),
        #                                      progressbar=True,
        #                                      chunk_size=10000))

    intint_reg = u.Quantity(intint_reg) * (u.solMass / u.pc**2)
    lwidth_reg = u.Quantity(lwidth_reg) * (u.m / u.s)

    intint_bins.append(intint_reg)
    lwidth_bins.append(lwidth_reg)

snr_props = {"shell_rad": [], "vol": [], "v_exp": [], "hole_mass": [],
             "shell_mass": [], "vol_dens": [], "energy": []}

# Half bins except for the first at 0.
dist_bin_corr = u.Quantity([dist_limit[0].value] + list(dist_limit[1:].value - 50)) * u.pc

show_plots = True

for i, (obj, obj2) in enumerate(zip(intint_bins, lwidth_bins)):

    out_props = find_bubble_props(dist_bin_corr, obj[:, 0], obj2[:, 0],
                                  snr_tab['D'][i] * u.pc)

    snr_props['shell_rad'].append(out_props[0].value)
    snr_props['vol'].append(out_props[1].value)
    snr_props['v_exp'].append(out_props[2].value)
    snr_props['hole_mass'].append(out_props[3].value)
    snr_props['shell_mass'].append(out_props[4].value)
    snr_props['vol_dens'].append(out_props[5].value)
    snr_props['energy'].append(out_props[6].value)

    if show_plots:
        fig = plt.figure(figsize=(12, 6))
        fig.add_subplot(131)
        pix = np.where(index_mask == i + 1)
        xlow, xhigh = np.min(pix[1]), np.max(pix[1])
        ylow, yhigh = np.min(pix[0]), np.max(pix[0])
        lim_slice = [slice(ylow - 50, yhigh + 50), slice(xlow - 50, xhigh + 50)]
        plt.imshow(mom0.value[lim_slice], origin='lower')
        plt.contour((index_mask == i + 1)[lim_slice], colors='b')
        fig.add_subplot(132)
        plt.errorbar(dist_limit.value, obj[:, 0].value, yerr=obj[:, 1].value, drawstyle='steps-mid')
        plt.xlabel("Distance (pc)")
        plt.ylabel(r"Surf. Density (Msol/pc$^2$)")
        fig.add_subplot(133)
        plt.errorbar(dist_limit.value, obj2[:, 0].value / 1000.,
                     yerr=obj2[:, 1].value / 1000., drawstyle='steps-mid')
        plt.xlabel("Distance (pc)")
        plt.ylabel(r"Line Width (km/s)")
        plt.tight_layout()
        plt.draw()
        print(out_props)
        input("{}".format(i + 1))
        plt.close()

snr_props['shell_rad'] = snr_props['shell_rad'] * u.pc
snr_props['vol'] = snr_props['vol'] * u.pc**3
snr_props['v_exp'] = snr_props['v_exp'] * u.km / u.s
snr_props['hole_mass'] = snr_props['hole_mass'] * u.solMass
snr_props['shell_mass'] = snr_props['shell_mass'] * u.solMass
snr_props['vol_dens'] = snr_props['vol_dens'] * u.cm**-3
snr_props['energy'] = snr_props['energy'] * u.erg

# Now we want to do something similar around the O-stars
# Using colour cuts for the half-brick near 604 (for now)

# import dask.dataframe as dd
import pandas as pd

df_phot = pd.read_hdf("/Volumes/Travel_Data/M33_2/Hubble/14610_M33-B01_1.phot.Ocut.hdf5",
                      key='data')

index_mask_O, dist_transf_O = \
    distance_map_from_catalogue(gal, df_phot, hi_cube.header,
                                diam_key=None, ra_key='ra',
                                dec_key='dec')

stacked_spectra_O = []
lwidth_bins_O = []
intint_bins_O = []

# Pick out individual regions
# labels, num = nd.label(dist_transf_O == 0.)
num = index_mask_O.max()

for n in ProgressBar(range(1, num + 1)):

    reg_mask = index_mask_O == n

    dist_transf_reg = nd.distance_transform_edt(~reg_mask)

    lwidth_reg = []
    intint_reg = []

    # Calculate avg properties within the region
    lwidth_reg.append([np.nanmean(lwidth[reg_mask].value),
                       np.nanstd(lwidth[reg_mask].value) / np.sqrt(reg_mask.sum() / 41.)])
    intint_reg.append([np.nanmean(moment0_surfdens[reg_mask].value),
                       np.nanstd(moment0_surfdens[reg_mask].value) / np.sqrt(reg_mask.sum() / 41.)])

    for i, (low, high) in enumerate(zip(dist_limit[:-1], dist_limit[1:])):

        # print("On bin {}".format(i + 1))

        dist_ang_low = (low / gal.distance.to(u.pc)).value * u.rad
        dist_pix_low = dist_ang_low.to(u.deg).value / np.abs(hi_cube.header['CDELT2'])

        dist_ang_high = (high / gal.distance.to(u.pc)).value * u.rad
        dist_pix_high = dist_ang_high.to(u.deg).value / np.abs(hi_cube.header['CDELT2'])

        dist_mask = np.logical_and(dist_transf_reg > dist_pix_low,
                                   dist_transf_reg <= dist_pix_high)

        num_beams = dist_mask.sum() / 41.

        intint_reg.append([np.nanmean(moment0_surfdens[dist_mask].value),
                           np.nanstd(moment0_surfdens[dist_mask].value) / np.sqrt(num_beams)])
        lwidth_reg.append([np.nanmean(lwidth[dist_mask].value),
                           np.nanstd(lwidth[dist_mask].value) / np.sqrt(num_beams)])

        # stacked_spectra.append(stack_spectra(hi_cube, peak_vel,
        #                                      xy_posns=np.where(dist_mask),
        #                                      progressbar=True,
        #                                      chunk_size=10000))

    intint_reg = u.Quantity(intint_reg) * (u.solMass / u.pc**2)
    lwidth_reg = u.Quantity(lwidth_reg) * (u.m / u.s)

    intint_bins_O.append(intint_reg)
    lwidth_bins_O.append(lwidth_reg)


# Plot individual profiles

o_props = {"shell_rad": [], "vol": [], "v_exp": [], "hole_mass": [],
             "shell_mass": [], "vol_dens": [], "energy": []}

# Half bins except for the first at 0.
dist_bin_corr = u.Quantity([dist_limit[0].value] + list(dist_limit[1:].value - 50)) * u.pc

show_plots = False

for i, (obj, obj2) in enumerate(zip(intint_bins, lwidth_bins)):

    out_props = find_bubble_props(dist_bin_corr, obj[:, 0], obj2[:, 0],
                                  0. * u.pc)

    o_props['shell_rad'].append(out_props[0].value)
    o_props['vol'].append(out_props[1].value)
    o_props['v_exp'].append(out_props[2].value)
    o_props['hole_mass'].append(out_props[3].value)
    o_props['shell_mass'].append(out_props[4].value)
    o_props['vol_dens'].append(out_props[5].value)
    o_props['energy'].append(out_props[6].value)

    if show_plots:
        plt.subplot(131)
        pix = np.where(index_mask_O == i + 1)
        if len(pix[0]) == 0:
            print("Found duplicated pixel location for {}.".format(i))
            continue
        xlow, xhigh = np.min(pix[1]), np.max(pix[1])
        ylow, yhigh = np.min(pix[0]), np.max(pix[0])
        lim_slice = [slice(ylow - 50, yhigh + 50), slice(xlow - 50, xhigh + 50)]
        plt.imshow(mom0.value[lim_slice], origin='lower')
        plt.contour((index_mask_O == i + 1)[lim_slice], colors='b')
        plt.subplot(132)
        plt.errorbar(dist_limit.value, obj[:, 0].value, yerr=obj[:, 1].value, drawstyle='steps-mid')
        plt.subplot(133)
        plt.errorbar(dist_limit.value, obj2[:, 0].value, yerr=obj2[:, 1].value, drawstyle='steps-mid')
        plt.draw()
        print(out_props)
        input("{}".format(i + 1))
        plt.clf()

o_props['shell_rad'] = o_props['shell_rad'] * u.pc
o_props['vol'] = o_props['vol'] * u.pc**3
o_props['v_exp'] = o_props['v_exp'] * u.km / u.s
o_props['hole_mass'] = o_props['hole_mass'] * u.solMass
o_props['shell_mass'] = o_props['shell_mass'] * u.solMass
o_props['vol_dens'] = o_props['vol_dens'] * u.cm**-3
o_props['energy'] = o_props['energy'] * u.erg

# Convert into dataframes and save
snr_hi_tab = Table([Column(snr_props[key]) for key in snr_props],
                names=snr_props.keys())
snr_hi_tab.write("/Volumes/Travel_Data/M33_2/HI/HI_snr_props.csv")

o_tab = Table([Column(o_props[key]) for key in o_props],
                names=o_props.keys())

o_tab.write("/Volumes/Travel_Data/M33_2/HI/HI_Ostar_props_M33-B01_1.csv")

# Compare some properties together

# Define the save path
import os
allfigs_path = lambda x: os.path.join(os.path.expanduser("~/Dropbox/Various Plots/M33/"), x)

_ = plt.hist(snr_hi_tab['shell_rad'], bins='auto', alpha=0.3, label='SNR')
_ = plt.hist(o_tab['shell_rad'], bins='auto', alpha=0.3, label='O')
plt.xlabel("Shell Radius (pc)")
plt.tight_layout()
plt.savefig(allfigs_path("feedback/feedback_sources_HI_shell_radius.png"))
plt.savefig(allfigs_path("feedback/feedback_sources_HI_shell_radius.pdf"))
plt.close()

_ = plt.hist(np.log10(snr_hi_tab['energy']
             [np.isfinite(snr_hi_tab['energy']) & (snr_hi_tab['energy'] > 0.)]),
             bins='auto', alpha=0.3, label='SNR')
_ = plt.hist(np.log10(o_tab['energy']
             [np.isfinite(o_tab['energy']) & (o_tab['energy'] > 0.)]),
             bins='auto', alpha=0.3, label='O')
plt.xlabel("log Energy (erg)")
plt.tight_layout()
plt.savefig(allfigs_path("feedback/feedback_sources_HI_energy.png"))
plt.savefig(allfigs_path("feedback/feedback_sources_HI_energy.pdf"))
plt.close()

# Properties  are similar. Can probably only believe a handful of small
# (<200 pc) but resolved sources.