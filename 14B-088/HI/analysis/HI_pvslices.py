
'''
Create PV-slices of the spectrally down-sampled HI cube.
'''

from spectral_cube import SpectralCube, Projection
import pvextractor as pv
from astropy.utils.console import ProgressBar
from astropy.io import fits
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt

from paths import fourteenB_HI_data_wGBT_path, fourteenB_wGBT_HI_file_dict
from galaxy_params import gal_feath as gal
from constants import distance


# Radius cut-off from galaxy frame to observation frame
def obs_radius(radius, PA, gal):

    ang_term = (np.cos(PA))**2 + (np.sin(PA) / np.cos(gal.inclination))**2

    return np.sqrt(radius**2 / ang_term)


def phys_to_ang(phys_size, distance=distance):
    '''
    Convert from angular to physical scales
    '''
    return (phys_size.to(distance.unit) / distance) * u.rad


if __name__ == "__main__":

    cube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.1kms.fits"))
    mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict["Moment0"])[0])

    pa_angles = gal.position_angles(header=mom0.header)
    radius = gal.radius(header=mom0.header)

    # Set the angles to create PV slices from. Angles defined wrt to the disk PA.
    # The cuts at  PA - 165 and 175 are to match the slices presented in Kam+17
    thetas = np.array([0, 45, 90, 135, gal.position_angle.value - 165,
                       gal.position_angle.value - 175]) * u.deg
    max_rad = 8.5 * u.kpc


    check_paths = False

    pvslices = []
    paths = []
    for theta in ProgressBar(thetas):

        # Adjust path length based on inclination
        obs_rad = obs_radius(max_rad, theta, gal)

        # Now convert the physical size in the observed frame to an angular size
        ang_length = 2 * phys_to_ang(obs_rad, distance)

        ang_width = 2 * phys_to_ang(obs_radius(max_rad, theta + 90 * u.deg, gal),
                                    distance)

        pv_path = pv.PathFromCenter(gal.center_position, length=ang_length,
                                    angle=theta + gal.position_angle,
                                    sample=20,
                                    width=ang_width)
        paths.append(pv_path)

        if check_paths:

            plt.imshow(mom0.value, origin='lower')
            plt.contour(radius <= max_rad, colors='r')
            center = gal.to_center_position_pixel(wcs=cube.wcs)
            plt.plot(center[0], center[1], 'bD')

            for i, posn in enumerate(pv_path.get_xy(cube.wcs)):
                if i == 0:
                    symb = "c^"
                else:
                    symb = "g^"
                plt.plot(posn[0], posn[1], symb)

            plt.draw()
            raw_input("{}".format(theta))
            plt.clf()

        # Set NaNs to zero. We're averaging over very large areas here.
        pvslice = pv.extract_pv_slice(cube, pv_path, respect_nan=False)

        filename = "downsamp_1kms/M33_14B-088_HI.GBT_feathered_PA_{}_pvslice.fits".format(int(theta.value))
        pvslice.writeto(fourteenB_HI_data_wGBT_path(filename, no_check=True), overwrite=True)

        pvslices.append(pvslice)

    # Make a major axis PV slice with the rotation subtracted cube
    rotcube = SpectralCube.read(fourteenB_HI_data_wGBT_path("downsamp_1kms/M33_14B-088_HI.clean.image.GBT_feathered.rotation_corrected.1kms.fits"))

    theta = 0.0 * u.deg
    obs_rad = obs_radius(max_rad, theta, gal)

    # Now convert the physical size in the observed frame to an angular size
    ang_length = 2 * phys_to_ang(obs_rad, distance)

    ang_width = 2 * phys_to_ang(obs_radius(max_rad, theta + 90 * u.deg, gal),
                                distance)

    pv_path = pv.PathFromCenter(gal.center_position, length=ang_length,
                                angle=theta + gal.position_angle,
                                sample=20,
                                width=ang_width)

    pvslice = pv.extract_pv_slice(rotcube, pv_path, respect_nan=False)

    filename = "downsamp_1kms/M33_14B-088_HI.GBT_feathered.rotation_corrected_PA_{}_pvslice.fits".format(int(theta.value))
    pvslice.writeto(fourteenB_HI_data_wGBT_path(filename, no_check=True), overwrite=True)


    # Make thinner PV slices with the normal cube
    pvslices_thin = []
    paths_thin = []
    for theta in ProgressBar(thetas):

        # Adjust path length based on inclination
        obs_rad = obs_radius(max_rad, theta, gal)

        # Now convert the physical size in the observed frame to an angular size
        ang_length = 2 * phys_to_ang(obs_rad, distance)

        ang_width = 200

        pv_path = pv.PathFromCenter(gal.center_position, length=ang_length,
                                    angle=theta + gal.position_angle,
                                    sample=20,
                                    width=ang_width)
        paths_thin.append(pv_path)

        if check_paths:

            plt.imshow(mom0.value, origin='lower')
            plt.contour(radius <= max_rad, colors='r')
            center = gal.to_center_position_pixel(wcs=cube.wcs)
            plt.plot(center[0], center[1], 'bD')

            for i, posn in enumerate(pv_path.get_xy(cube.wcs)):
                if i == 0:
                    symb = "c^"
                else:
                    symb = "g^"
                plt.plot(posn[0], posn[1], symb)

            plt.draw()
            raw_input("{}".format(theta))
            plt.clf()

        # Set NaNs to zero. We're averaging over very large areas here.
        pvslice = pv.extract_pv_slice(cube, pv_path, respect_nan=False)

        filename = "downsamp_1kms/M33_14B-088_HI.GBT_feathered_PA_{}_pvslice_200arcsec_width.fits".format(int(theta.value))
        pvslice.writeto(fourteenB_HI_data_wGBT_path(filename, no_check=True), overwrite=True)

        pvslices_thin.append(pvslice)
