
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

def generate_vrot_model():
    data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"


    data = \
        Table.read(os.path.join(data_path,
                                'diskfit_noasymm_nowarp_output/rad.out.csv'))

    pars, pcov = curve_fit(vcirc, data['r'], data['Vt'], sigma=data['eVt'],
                           absolute_sigma=True, p0=(1., 100., 1000.))

    print("n: {0} +/- {1}".format(pars[0], np.sqrt(pcov[0, 0])))
    print("vmax: {0} +/- {1}".format(pars[1], np.sqrt(pcov[1, 1])))
    print("rmax: {0} +/- {1}".format(pars[2], np.sqrt(pcov[2, 2])))

    return pars, pcov

def return_smooth_model(header):

    gal = Galaxy("M33")

    radii = gal.radius(header=header).value
    pas = gal.position_angles(header=header).value

    scale = wcs.utils.proj_plane_pixel_scales(wcs.WCS(header))[0]
    # Distance scaling (1" ~ 4 pc). Conversion is deg to kpc
    dist_scale = (np.pi / 180.) * gal.distance.to(u.pc).value

    pars, pcov = generate_vrot_model()

    # Convert rmax to pc
    mod_pars = pars.copy()
    mod_pars[2] *= scale * dist_scale

    # Put into m/s.
    smooth_model = (vcirc(radii, *mod_pars) * 1000.) * np.cos(pas) * \
        np.sin(gal.inclination).value

    # Shift by Vsys (m / s)
    vsys = -180610.
    smooth_model += vsys

    return smooth_model

if __name__ == "__main__":
    make_plot = True
    make_rotmodel = False


    pars, pcov = generate_vrot_model()

    data_path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"

    mom1_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom1.fits"
    # mom1_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.3.ellip_mask.mom1.fits"
    mom1 = \
        fits.open(os.path.join(data_path, mom1_name))

    # Pixel scales (1 ~ 3")
    scale = wcs.utils.proj_plane_pixel_scales(wcs.WCS(mom1[0].header))[0]
    # Distance scaling (1" ~ 4 pc). Conversion is deg to kpc
    dist_scale = (np.pi / 180.) * 840.

    if make_plot:

        data = \
            Table.read(os.path.join(data_path,
                                    'diskfit_noasymm_nowarp_output/rad.out.csv'))

        phys_radius = data["r"] * scale * dist_scale
        plot_pars = pars.copy()
        plot_pars[2] *= scale * dist_scale
        p.errorbar(phys_radius, data["Vt"], yerr=data['eVt'],
                   fmt='-', color='b', label="This work", drawstyle='steps-mid')
        p.plot(phys_radius, vcirc(phys_radius, *plot_pars), 'r-')
        p.ylabel(r"V$_{\mathrm{circ}}$ (km / s)")
        p.xlabel(r"Radius (kpc)")
        p.grid()
        p.draw()

        raw_input("Next plot?")
        p.clf()

        # load in the Corbelli curve for comparison
        corbelli = Table.read(os.path.expanduser("~/Dropbox/code_development/VLA_Lband/14B-088/HI/analysis/rotation_curves/corbelli_rotation_curve.csv"))

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

        p.ion()

    if make_rotmodel:

        smooth_model = return_smooth_model(mom1[0].header)

        # Save the smooth model.
        new_header = mom1[0].header.copy()
        new_header["COMMENT"] = "Smooth rotation model of DISKFIT output. " \
            "Uses Eq.5 from Meidt et al. 2008. n={0:.2f}+/-{1:.2f}, " \
            "Vmax={2:.2f}+/-{3:.2f} km/s, rmax={4:.2f}+/-{5:.2f} pix".\
            format(pars[0], np.sqrt(pcov[0, 0]), pars[1], np.sqrt(pcov[1, 1]),
                   pars[2], np.sqrt(pcov[2, 2]))
        new_header["BUNIT"] = "m / s"
        new_hdu = fits.PrimaryHDU(smooth_model, header=new_header)
        new_hdu.writeto(os.path.join(data_path,
                                     "diskfit_noasymm_nowarp_output/rad.fitmod.fits"),
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
        new_hdu.writeto(os.path.join(data_path,
                                     "diskfit_noasymm_nowarp_output/rad.fitres.fits"),
                        clobber=True)
