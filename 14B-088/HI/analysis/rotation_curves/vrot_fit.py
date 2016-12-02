
import numpy as np
from scipy.optimize import curve_fit
from astropy.table import Table
import os
from astropy.io import fits
import astropy.units as u
import astropy.wcs as wcs
import matplotlib.pyplot as p
from galaxies import Galaxy

'''
Fit Meidt+08 Eq. 5 (Eq. 1 Faber & Gallagher 79)
'''


def vcirc(r, *pars):
    n, vmax, rmax = pars
    numer = vmax * (r / rmax)
    denom = np.power((1 / 3.) + (2 / 3.) *
                     np.power(r / rmax, n), (3 / (2 * n)))
    return numer / denom


def generate_vrot_model(table_name, verbose=True):
    '''
    Parameters
    ----------
    table_name : str
        Name and path of the csv table produced by `run_diskfit.py`
    '''

    if isinstance(table_name, basestring):
        data = Table.read(table_name)
    else:
        data = table_name

    pars, pcov = curve_fit(vcirc, data['r'], data['Vt'], sigma=data['eVt'],
                           absolute_sigma=True, p0=(1., 100., 1000.))

    if verbose:
        print("n: {0} +/- {1}".format(pars[0], np.sqrt(pcov[0, 0])))
        print("vmax: {0} +/- {1}".format(pars[1], np.sqrt(pcov[1, 1])))
        print("rmax: {0} +/- {1}".format(pars[2], np.sqrt(pcov[2, 2])))

    return pars, pcov


def return_smooth_model(table_name, header, gal):

    radii = gal.radius(header=header).value
    pas = gal.position_angles(header=header).value

    scale = wcs.utils.proj_plane_pixel_scales(wcs.WCS(header))[0]
    # Distance scaling (1" ~ 4 pc). Conversion is deg to kpc
    dist_scale = (np.pi / 180.) * gal.distance.to(u.pc).value

    pars, pcov = generate_vrot_model(table_name)

    # Convert rmax to pc
    mod_pars = pars.copy()
    mod_pars[2] *= scale * dist_scale

    # Put into m/s.
    smooth_model = (vcirc(radii, *mod_pars) * 1000.) * np.cos(pas) * \
        np.sin(gal.inclination).value

    # Shift by Vsys (m / s)
    smooth_model += gal.vsys.to(u.m / u.s).value

    return smooth_model


def update_galaxy_params(gal, param_table, spatial_wcs):
    '''
    Use the fit values from fit rather than the hard-coded values in galaxies.
    '''

    from astropy.coordinates import Angle, SkyCoord

    gal.inclination = Angle(param_table["inc"] * u.deg)[0]
    gal.position_angle = Angle(param_table["PA"] * u.deg)[0]
    gal.vsys = (param_table["Vsys"] * u.km / u.s)[0]

    # The positions in the table are in pixels, so convert to the sky using
    # the spatial WCS info.
    xcent, ycent = param_table["xcent"], param_table["ycent"]
    ra_cent, dec_cent = spatial_wcs.celestial.wcs_pix2world(xcent, ycent, 0)

    gal.center_position = SkyCoord(ra_cent, dec_cent, unit=(u.deg, u.deg),
                                   frame='fk5')


if __name__ == "__main__":

    import sys

    from analysis.paths import (fourteenB_HI_data_path, paper1_figures_path,
                                c_hi_analysispath)

    make_plot = True
    make_rotmodel = True

    gal = Galaxy("M33")

    try:
        folder_name = sys.argv[1]
    except IndexError:

        diskfit_runs = [dir for dir in os.listdir(fourteenB_HI_data_path(""))
                        if dir.startswith("diskfit")]

        print("Available diskfit runs: " + str(diskfit_runs))

        folder_name = raw_input("Give folder name of the diskfit run: ")

    # My convention has been diskfit_params_output. Only keep the params parts
    # for the names of the output plots.
    params = "_".join(folder_name.split("_")[1:-1])

    # Load in moment 1
    # mom1_name = fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3.ellip_mask.mom1.fits")
    mom1_name = fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3.ellip_mask.mom1.fits")
    mom1 = fits.open(mom1_name)
    mom1_wcs = wcs.WCS(mom1[0].header)

    # Since diskfit gives updated galaxy parameters, read them in and use
    # instead of the hard-wired values in Galaxy
    param_name = \
        fourteenB_HI_data_path("{}/rad.out.params.csv".format(folder_name))

    param_table = Table.read(param_name)

    update_galaxy_params(gal, param_table, mom1_wcs)

    table_name = \
        fourteenB_HI_data_path("{}/rad.out.csv".format(folder_name))

    data = Table.read(table_name)

    pars, pcov = generate_vrot_model(table_name)

    # Pixel scales (1 ~ 3")
    scale = wcs.utils.proj_plane_pixel_scales(mom1_wcs)[0]
    # Distance scaling (1" ~ 4 pc). Conversion is deg to kpc
    dist_scale = (np.pi / 180.) * gal.distance.to(u.kpc).value

    if make_plot:

        phys_radius = data["r"] * scale * dist_scale
        plot_pars = pars.copy()
        plot_pars[2] *= scale * dist_scale
        p.errorbar(phys_radius, data["Vt"], yerr=data['eVt'],
                   fmt='-', color='b', label="This work",
                   drawstyle='steps-mid')
        p.plot(phys_radius, vcirc(phys_radius, *plot_pars), 'r-')
        p.ylabel(r"V$_{\mathrm{circ}}$ (km / s)")
        p.xlabel(r"Radius (kpc)")
        p.grid()
        p.draw()
        p.savefig(paper1_figures_path("M33_vrot_{}_wfit.pdf".format(params)))
        p.savefig(paper1_figures_path("M33_vrot_{}_wfit.png".format(params)))

        # raw_input("Next plot?")
        p.clf()

        # load in the Corbelli curve for comparison
        corbelli = Table.read(c_hi_analysispath("rotation_curves/corbelli_rotation_curve.csv"))

        p.errorbar(phys_radius, data["Vt"], yerr=data['eVt'],
                   fmt='-', color='b', label="This work", drawstyle='steps-mid')
        p.errorbar(corbelli["R"][corbelli["R"] <= 8.0],
                   corbelli["Vr"][corbelli["R"] <= 8.0],
                   yerr=corbelli["dVr"][corbelli["R"] <= 8.0],
                   fmt='--', color='r', label="Corbelli et al. (2014)",
                   drawstyle='steps-mid')
        p.ylabel(r"V$_{\mathrm{circ}}$ (km / s)")
        p.xlabel(r"Radius (kpc)")
        p.legend(loc='lower right')
        p.grid()
        p.draw()

        p.savefig(paper1_figures_path("M33_vrot_{}_wCorbelli.pdf".format(params)))
        p.savefig(paper1_figures_path("M33_vrot_{}_wCorbelli.png".format(params)))

        p.close()
        p.ion()

    if make_rotmodel:

        # First make version with the ellipse masked moment 1 loaded in above.
        smooth_model = return_smooth_model(data, mom1[0].header, gal)

        # Save the smooth model.
        new_header = mom1[0].header.copy()
        new_header["COMMENT"] = "Smooth rotation model of DISKFIT output. " \
            "Uses Eq.5 from Meidt et al. 2008. n={0:.2f}+/-{1:.2f}, " \
            "Vmax={2:.2f}+/-{3:.2f} km/s, rmax={4:.2f}+/-{5:.2f} pix".\
            format(pars[0], np.sqrt(pcov[0, 0]), pars[1], np.sqrt(pcov[1, 1]),
                   pars[2], np.sqrt(pcov[2, 2]))
        new_header["BUNIT"] = "m / s"
        new_hdu = fits.PrimaryHDU(smooth_model, header=new_header)
        new_hdu.writeto(fourteenB_HI_data_path("{}/rad.fitmod_ellipmask.fits".format(folder_name),
                                               no_check=True),
                        clobber=True)

        # And the non-circular residuals
        new_header = mom1[0].header.copy()
        new_header["COMMENT"] = "Residuals from rotation model of DISKFIT output. " \
            "Uses Eq.5 from Meidt et al. 2008. n={0:.2f}+/-{1:.2f}, " \
            "Vmax={2:.2f}+/-{3:.2f} km/s, rmax={4:.2f}+/-{5:.2f} pix".\
            format(pars[0], np.sqrt(pcov[0, 0]), pars[1], np.sqrt(pcov[1, 1]),
                   pars[2], np.sqrt(pcov[2, 2]))
        new_header["BUNIT"] = "m / s"
        new_hdu = fits.PrimaryHDU(mom1[0].data - smooth_model, header=new_header)
        new_hdu.writeto(fourteenB_HI_data_path("{}/rad.fitres_ellipmask.fits".format(folder_name),
                                               no_check=True),
                        clobber=True)

        # And now with the normal pbcov masked grid.
        mom1_name = fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom1.fits")
        mom1 = fits.open(mom1_name)

        smooth_model = return_smooth_model(data, mom1[0].header, gal)

        # Save the smooth model.
        new_header = mom1[0].header.copy()
        new_header["COMMENT"] = "Smooth rotation model of DISKFIT output. " \
            "Uses Eq.5 from Meidt et al. 2008. n={0:.2f}+/-{1:.2f}, " \
            "Vmax={2:.2f}+/-{3:.2f} km/s, rmax={4:.2f}+/-{5:.2f} pix".\
            format(pars[0], np.sqrt(pcov[0, 0]), pars[1], np.sqrt(pcov[1, 1]),
                   pars[2], np.sqrt(pcov[2, 2]))
        new_header["BUNIT"] = "m / s"
        new_hdu = fits.PrimaryHDU(smooth_model, header=new_header)
        new_hdu.writeto(fourteenB_HI_data_path("{}/rad.fitmod.fits".format(folder_name),
                                               no_check=True),
                        clobber=True)

        # And the non-circular residuals
        new_header = mom1[0].header.copy()
        new_header["COMMENT"] = "Residuals from rotation model of DISKFIT output. " \
            "Uses Eq.5 from Meidt et al. 2008. n={0:.2f}+/-{1:.2f}, " \
            "Vmax={2:.2f}+/-{3:.2f} km/s, rmax={4:.2f}+/-{5:.2f} pix".\
            format(pars[0], np.sqrt(pcov[0, 0]), pars[1], np.sqrt(pcov[1, 1]),
                   pars[2], np.sqrt(pcov[2, 2]))
        new_header["BUNIT"] = "m / s"
        new_hdu = fits.PrimaryHDU(mom1[0].data - smooth_model, header=new_header)
        new_hdu.writeto(fourteenB_HI_data_path("{}/rad.fitres.fits".format(folder_name),
                                               no_check=True),
                        clobber=True)
