
from astropy.io import fits
from spectral_cube import SpectralCube, VaryingResolutionSpectralCube
from spectral_cube.lower_dimensional_structures import OneDSpectrum
import astropy.units as u
import numpy as np
import os
from astropy.utils.console import ProgressBar
from itertools import izip, repeat
from turbustat.statistics.stats_utils import fourier_shift


def spectrum_shifter(spectrum, v0, vcent):
    '''
    Shift the central velocity of a spectrum by the difference if v0 and vcent.

    Parameters
    ----------
    spectrum : `~spectral_cube.lower_dimensional_objects.OneDSpectrum`
        1D spectrum to shift.
    v0 : `~astropy.units.Quantity`
        Velocity to shift spectrum to.
    vcent : `~astropy.units.Quantity`
        Initial center velocity of the spectrum.
    '''

    vdiff = np.abs(spectrum.spectral_axis[1] - spectrum.spectral_axis[0])
    vel_unit = vdiff.unit

    pix_shift = (vcent.to(vel_unit) - v0.to(vel_unit)) / vdiff

    shifted = fourier_shift(spectrum, pix_shift)

    if hasattr(spectrum, "beams"):
        beams = spectrum.beams
    else:
        beams = None

    return OneDSpectrum(shifted, unit=spectrum.unit, wcs=spectrum.wcs,
                        meta=spectrum.meta, spectral_unit=vel_unit,
                        beams=beams)


def mulproc_spectrum_shifter(inputs):
    y, x, spec, vcent, v0 = inputs

    return spectrum_shifter(spec, v0, vcent), y, x


def cube_shifter(cube, velocity_surface, v0=None, save_shifted=False,
                 save_name=None, xy_posns=None, pool=None,
                 return_spectra=True, chunk_size=20000, is_mask=False,
                 verbose=False):
    '''
    Shift spectra in a cube according to a given velocity surface (peak
    velocity, centroid, rotation model, etc.).
    '''

    if not save_shifted and not return_spectra:
        raise Exception("One of 'save_shifted' or 'return_spectra' must be "
                        "enabled.")

    if not np.isfinite(velocity_surface).any():
        raise Exception("velocity_surface contains no finite values.")

    if xy_posns is None:
        # Only compute where a shift can be found
        xy_posns = np.where(np.isfinite(velocity_surface))

    if v0 is None:
        # Set to near the center velocity of the cube if not given.
        v0 = cube.spectral_axis[cube.shape[0] // 2]
    else:
        if not isinstance(v0, u.Quantity):
            raise u.UnitsError("v0 must be a quantity.")
        spec_unit = cube.spectral_axis.unit
        if not v0.unit.is_equivalent(spec_unit):
            raise u.UnitsError("v0 must have units equivalent to the cube's"
                               " spectral unit ().".format(spec_unit))

    # Adjust the header to have velocities centered at v0.
    new_header = cube.header.copy()
    new_header["CRVAL3"] = new_header["CRVAL3"] - v0.to(u.m / u.s).value

    if save_shifted:

        if os.path.exists(save_name):
            raise Exception("The file name {} already "
                            "exists".format(save_name))

        output_fits = fits.StreamingHDU(save_name, new_header)

        if is_mask:
            fill_plane = np.zeros_like(cube[0].value, dtype=np.int16)
        else:
            fill_plane = np.zeros_like(cube[0].value) * np.NaN

        for chan in xrange(cube.shape[0]):
            output_fits.write(fill_plane)
        output_fits.close()

        output_fits = fits.open(save_name, mode='update')

    if return_spectra:
        all_shifted_spectra = []
        out_posns = []

    # Create chunks of spectra for read-out.
    nchunks = 1 + xy_posns[0].size // chunk_size

    if verbose:
        iterat = ProgressBar(xrange(nchunks))
    else:
        iterat = xrange(nchunks)

    for i in iterat:

        y_posns = xy_posns[0][i * chunk_size:(i + 1) * chunk_size]
        x_posns = xy_posns[1][i * chunk_size:(i + 1) * chunk_size]

        gen = ((y, x, cube[:, y, x], velocity_surface[y, x], v0) for y, x in
               izip(y_posns, x_posns))

        if pool is not None:
            shifted_spectra = pool.map(mulproc_spectrum_shifter, gen)
        else:
            shifted_spectra = map(mulproc_spectrum_shifter, gen)

        if save_shifted:
            for out in shifted_spectra:
                if is_mask:
                    spec = (out[0].value > 0.5).astype(np.int)
                else:
                    spec = out[0].value
                output_fits[0].data[:, out[1], out[2]] = spec

            output_fits.flush()

        if return_spectra:
            all_shifted_spectra.extend([out[0] for out in shifted_spectra])
            out_posns.extend([out[1:] for out in shifted_spectra])

    if save_shifted:

        output_fits.flush()
        output_fits.close()

        # Append the beam table onto the output file.
        if isinstance(cube, VaryingResolutionSpectralCube):
            from spectral_cube.cube_utils import beams_to_bintable
            output_fits = fits.open(save_name, mode='append')
            output_fits.append(beams_to_bintable(cube.beams))
            output_fits.flush()
            output_fits.close()

    if return_spectra:
        return all_shifted_spectra, out_posns


def test_shifter(shape=(100, 100, 100), sigma=8., amp=1.):
    '''
    Use a set of identical Gaussian profiles randomly offset to ensure the
    shifted spectrum has the correct properties.
    '''

    def gaussian(x, amp, mean, sigma):
        return amp * np.exp(- (x - mean)**2 / (2 * sigma**2))

    test_cube = np.empty(shape)
    mean_positions = np.empty(shape[1:])

    spec_middle = shape[0] / 2
    spec_quarter = shape[0] / 4

    np.random.seed(247825498)

    spec_inds = np.mgrid[-spec_middle:spec_middle]
    spat_inds = np.indices(shape[1:])
    for y, x in zip(spat_inds[0].flatten(), spat_inds[1].flatten()):

        # Lock the mean to within 25% from the centre
        mean_pos = np.random.uniform(low=-spec_quarter,
                                     high=spec_quarter)

        mean_positions[y, x] = mean_pos + spec_middle
        test_cube[:, y, x] = gaussian(spec_inds, amp, mean_pos, sigma)

    # Now just use fourier shift to shift all back to the centre
    test_shifted_cube = np.empty(shape)

    for y, x in zip(spat_inds[0].flatten(), spat_inds[1].flatten()):

        pixel_diff = spec_middle - mean_positions[y, x]
        test_shifted_cube[:, y, x] = \
            fourier_shift(test_cube[:, y, x], pixel_diff)

    # Now fit a Gaussian to the mean stacked profile.
    from scipy.optimize import curve_fit

    mean_stacked_profile = test_shifted_cube.mean(axis=(1, 2))

    fit_vals = curve_fit(gaussian, spec_inds, mean_stacked_profile)[0]

    np.testing.assert_allclose(fit_vals, np.array([amp, 0.0, sigma]),
                               atol=1e-3)
