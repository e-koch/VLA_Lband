
'''
Run a variety of turbulent statistics on the full M33 cube.

This will use a TON of memory. Recommend running on a cluster.

File structure setup to work on cedar in scratch space.

Saves the outputs for later use
'''

from spectral_cube import SpectralCube
from astropy.io import fits
from os.path import join as osjoin
import astropy.units as u
from astropy import log
import sys

from turbustat.statistics import (PowerSpectrum, VCA, VCS, PCA, SCF,
                                  StatMoments, DeltaVariance)


ncore = int(sys.argv[-1])

run_pspec = False
run_delvar = False
run_moments = False
run_vca = True
run_vcs = True
run_pca = True
run_scf = True

scratch_path = "/home/ekoch/scratch/M33_turbulence/"

cube_name = osjoin(scratch_path, "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.fits")
mom0_name = osjoin(scratch_path, "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.mom0.fits")

out_path = lambda x: osjoin(scratch_path, x)

if run_pspec:

    log.info("Running PowerSpectrum")

    mom0 = fits.open(mom0_name)[0]

    pspec = PowerSpectrum(mom0).run(use_pyfftw=True, threads=ncore - 1)

    pspec.save_results(out_path("pspec_m33_14B088.pkl"), keep_data=False)

    del pspec

if run_delvar:

    log.info("Running Delta-Variance")

    mom0 = fits.open(mom0_name)[0]

    delvar = DeltaVariance(mom0).run(use_pyfftw=True, threads=ncore - 1)

    delvar.save_results(out_path("delvar_m33_14B088.pkl"), keep_data=False)

    del delvar

if run_moments:

    log.info("Running Stat Moments")

    mom0 = fits.open(mom0_name)[0]

    # Run with a couple of different radii
    # Beam size (close to it)
    radii = [7, 14, 21] * u.pix

    for radius in radii:

        log.info("StatMoments radius: {}".format(radius))

        moments = StatMoments(mom0, radius=radius).run()

        moments.save_results(out_path("moments_m33_14B088_radius_{}pix.pkl".format(radius.value)),
                             keep_data=False)

        del moments

if run_vca:

    log.info("Running VCA")

    # Avoid loading in multiple times
    try:
        cube.shape
    except NameError:
        cube = SpectralCube.read(cube_name, memmap=False)

    # Run VCA over a number of channel sizes.
    # First one keeps the original channel width.
    chan_widths = [None, 0.4 * u.km / u.s, 0.6 * u.km / u.s,
                   0.8 * u.km / u.s, 1.0 * u.km / u.s,
                   1.6 * u.km / u.s, 2.0 * u.km / u.s,
                   3.0 * u.km / u.s, 4.0 * u.km / u.s,
                   5.0 * u.km / u.s, 6.0 * u.km / u.s,
                   8.0 * u.km / u.s, 10.0 * u.km / u.s,
                   20.0 * u.km / u.s, 40.0 * u.km / u.s]

    for chan in chan_widths:

        log.info("On VCA channel width {}".format(chan))

        vca = VCA(cube, channel_width=chan).run(use_pyfftw=True,
                                                threads=ncore - 1)

        if chan is None:
            chan = 0.2 * u.km / u.s

        vca.save_results(out_path("vca_m33_14B088_chanwidth_{}_kms.pkl".format(chan.value)),
                         keep_data=False)

        del vca


if run_vcs:

    log.info("Running VCS")

    # Avoid loading in multiple times
    try:
        cube.shape
    except NameError:
        cube = SpectralCube.read(cube_name, memmap=False)

    # Run VCS when varying the spatial resolution.
    # First one keeps the original beam size.

    from radio_beam import Beam

    majors = [None, 38 * u.arcsec, 57 * u.arcsec,
              76 * u.arcsec, 95 * u.arcsec]

    for major in majors:

        log.info("On VCS resolution {}".format(major))

        if major is None:
            conv_cube = cube
        else:
            new_beam = Beam(major)
            conv_cube = cube.convolve_to(new_beam)

        vcs = VCS(conv_cube, channel_width=chan).run(use_pyfftw=True,
                                                     threads=ncore - 1)

        if major is None:
            major = 19 * u.arcsec

        vcs.save_results(out_path("vcs_m33_14B088_chanwidth_{}_arcsec.pkl".format(major.value)),
                         keep_data=False)

        del vcs
        del conv_cube


if run_pca:

    log.info("Running PCA")

    # Avoid loading in multiple times
    try:
        cube.shape
    except NameError:
        cube = SpectralCube.read(cube_name, memmap=False)

    pca = PCA(cube, distance=840 * u.kpc).run(min_eigval=0.001,
                                              beam_fwhm=19 * u.arcsec,
                                              spatial_output_unit=u.pc,
                                              spectral_output_unit=u.km / u.s)

    pca.save_results(out_path("pca_m33_14B088.pkl"),
                     keep_data=False)

    del pca

if run_scf:

    log.info("Running SCF")

    # Note that the rolls are going to take a LONG time!!
    scf = SCF(cube, size=36).run()

    scf.save_results(out_path("scf_m33_14B088.pkl"), keep_data=False)

    del scf
