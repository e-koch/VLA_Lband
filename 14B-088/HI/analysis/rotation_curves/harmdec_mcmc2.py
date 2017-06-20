
'''
This is Dario Colombo's code.
'''

import warnings
import os
import emcee

import numpy as np
import math
import aplpy
# import fish

# from matplotlib import pyplot as plt
import matplotlib.pyplot as plt

from astropy import constants as consts
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy.table.table import Column
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.utils.console import ProgressBar


def clean_hd2d(hd):

    hd.update(NAXIS=2)

    for key in hd:
        if key[-1] == "3" or key[-1] == "4":
            hd.remove(key)

    return hd


def log_hd(p, hds, vobs, evobs, psi, incl, order):

    vmod = p[0]
    # priors =  ss.norm.logpdf(p[0],loc=hds[0],scale=20)

    for i in range(1, order + 1):
        vmod = vmod + p[2 * i - 1] * np.cos(i * psi) * np.sin(incl) \
            + p[2 * i] * np.sin(i * psi) * np.sin(incl)
        """
          priors = priors + ss.norm.logpdf(p[2*i-1],loc=hds[2*i-1],scale=20.) + \
            ss.norm.logpdf(p[2*i],loc=hds[2*i],scale=20.)
        """

    p1 = (vmod - vobs)**2 / evobs**2
    p1 = np.nansum(p1)

    lp = - p1  # + priors

    if np.isnan(lp):
        return -np.inf

    return lp


def deproj(data, hd, incl=45., pa=180., dist=5e6, aoff=None, doff=None,
           vsys=None):

    if aoff is None:
        aoff = hd.get('CRVAL1')

    if doff is None:
        doff = hd.get('CRVAL2')

    if vsys is None:
        vsys = np.nanmedian(data)

    w = wcs.WCS(hd)
    xint, yint = np.meshgrid(np.arange(data.shape[1]),
                             np.arange(data.shape[0]),
                             indexing='ij')

    ra, dec = w.all_pix2world(xint, yint, 0)

    # Physical properties and deprojection
    Distance = dist * u.pc
    PC2Deg = 1 / (np.pi * Distance / 180)
    CenterPosition = SkyCoord(aoff, doff, unit=(u.deg, u.deg), frame='fk5')
    PositionAngle = pa * u.deg
    Inclination = incl * u.deg
    Vsys = vsys * u.km / u.s
    Offsets = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    PAs = CenterPosition.position_angle(Offsets)
    GalPA = PAs - PositionAngle
    GCDist = Offsets.separation(CenterPosition)
    Rplane = Distance * np.tan(GCDist)
    Xplane = Rplane * np.cos(GalPA)
    Yplane = Rplane * np.sin(GalPA)
    Xgal = Xplane
    Ygal = Yplane / np.cos(Inclination)
    Rgal = (Xgal**2 + Ygal**2)**0.5
    Rgal = Rgal * PC2Deg * 3600
    Psigal = np.arctan2(Ygal, Xgal)

    return Rgal, Psigal


def harmdec(data, hd, props, order=1, blind=False):

    # Reading file and header
    data = data.squeeze()
    hd = clean_hd2d(hd)

    # Reading the properties
    # prop = [aoff, doff, dist, incl, pa, vsys]

    aoff = props[0]
    doff = props[1]
    dist = props[2]
    incl = props[3]
    pa = props[4]
    vsys = props[5]
    rmax = props[6]
    beam = props[7]

    # Deproject
    Rgal, Psigal = deproj(data, hd, incl=incl, pa=pa, dist=dist,
                          aoff=aoff, doff=doff, vsys=vsys)

    # Linearize the images
    rgal = Rgal.value.T
    psigal = Psigal.value.T
    rr = rgal.reshape(rgal.shape[0] * rgal.shape[1])
    ppsi = psigal.reshape(psigal.shape[0] * psigal.shape[1])
    vobs = data.reshape(data.shape[0] * data.shape[1])
    evobs = np.ones(vobs.shape)

    # Masking bad pixels
    rr[np.isnan(vobs)] = np.nan
    ppsi[np.isnan(vobs)] = np.nan

    # Kinematic parameters in radians
    incl = np.radians(incl)
    pa = np.radians(pa)

    # Parameters for the fit
    dr = beam
    rint = beam + beam / 2
    rmin = rint
    rmaj = rmin + dr

    if rmax is None:
        rmax = np.nanmax(rr)

    # Empty list to accomodate harmonic terms
    rads = []
    shds = [[] for i in range(2 * order + 1)]
    hds = [[] for i in range(2 * order + 1)]
    emhds = [[] for i in range(2 * order + 1)]
    ephds = [[] for i in range(2 * order + 1)]

    # Start the harmonic decomposition
    j = 0
    jmax = np.int((rmax - rmin) / beam)
    mask = np.zeros(vobs.shape)
    vmod = np.zeros(vobs.shape)
    rotmod = np.zeros(vobs.shape)

    # peixe = fish.ProgressFish(total=jmax)
    peixe = ProgressBar(jmax)

    if not blind:

            # plt.switch_backend('Qt4Agg')
            # plt.switch_backend('MacOSX')
        fig = plt.figure(figsize=(14, 6))
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)

    """
    ffig = plt.figure(figsize=(20,10))
    axf0 = ffig.add_subplot(131)
    axf1 = ffig.add_subplot(132)
    axf2 = ffig.add_subplot(133)
    """

    while rmaj < rmax:

        # peixe.animate(amount=j)
        peixe.update(j)

        # Select the ring
        idx = (rr > rmin) & (rr <= rmaj)
        rads.append(np.nanmedian(rr[idx]))
        svobs = vobs[idx]
        sevobs = evobs[idx]
        sppsi = ppsi[idx]

        mask[idx] = j

        # Guessing the harmonic terms with SVD
        A = np.zeros((len(svobs), 2 * order))
        for i in range(order):
            A[:, 2 * i] = np.cos((i + 1) * sppsi) * np.sin(incl)
            A[:, 2 * i + 1] = np.sin((i + 1) * sppsi) * np.sin(incl)

        B = svobs - vsys

        try:

            x, _, _, _ = np.linalg.lstsq(A, B)

            """
            U,s,V = np.linalg.svd(A, full_matrices = False)
            c = np.dot(U.T,B) # c = U^t*b
            w = np.linalg.solve(np.diag(s),c) # w = V^t*c
            x = np.dot(V.T,w)
            """

            shds[0].append(vsys)
            # emhds[0].append(1.)
            # ephds[0].append(1.)

            for i in range(1, order + 1):
                shds[2 * i - 1].append(x[2 * i - 2])
                shds[2 * i].append(x[2 * i - 1])
                # emhds[2*i-1].append(1.)
                # ephds[2*i-1].append(1.)
                # emhds[2*i].append(1.)
                # ephds[2*i].append(1.)

            # Optimize the model with MCMC
            ndim = 2 * order + 1
            nwalkers = 2 * ndim
            p0 = np.zeros((nwalkers, ndim))
            p0[:, 0] = np.random.randn(nwalkers) * 20 + vsys

            mhds = [0.]
            for i in range(1, ndim):
                p0[:, i] = np.random.randn(nwalkers) * 50 + shds[i][j]
                #p0[:,i] = np.random.uniform(size = nwalkers, low = shds[i][j]-20., high = shds[i][j]+20.)
                mhds.append(shds[i][j])

            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_hd,
                                            args=[mhds, svobs, sevobs, sppsi,
                                                  incl, order],
                                            threads=1)

            pos, prob, state = sampler.run_mcmc(p0, 500, thin=5)
            sampler.reset()
            pos, prob, state = sampler.run_mcmc(pos, 4000, thin=5)

            """
            axf0.plot(sampler.flatchain[:,0])
            axf1.plot(sampler.flatchain[:,1])
            axf2.plot(sampler.flatchain[:,2])
            plt.draw()
            """

            # Filling the hd parameter values
            for i in range(ndim):
                hds[i].append(shds[i][j])
                ephds[i].append(np.percentile(sampler.flatchain[:, i], 75) -
                                np.median(sampler.flatchain[:, i]))
                emhds[i].append(np.median(sampler.flatchain[:, i]) -
                                np.percentile(sampler.flatchain[:, i], 25))

        except ValueError:

            ndim = 2 * order + 1

            for i in range(ndim):
                shds[i].append(np.nan)
                hds[i].append(np.nan)
                ephds[i].append(np.nan)
                emhds[i].append(np.nan)

        # Making the model
        vmod[idx] = vsys
        for i in range(1, order + 1):
            vmod[idx] = vmod[idx] + hds[2 * i - 1][j] * np.cos(i * sppsi) * \
                np.sin(incl) + hds[2 * i][j] * np.sin(i * sppsi) * np.sin(incl)

        rotmod[idx] = hds[1][j] * np.cos(sppsi) * np.sin(incl)

        if not blind:

            vmodp = vmod.copy()
            vmodp[vmodp == 0] = np.nan
            vmodp = vmodp.reshape(data.shape[0], data.shape[1])
            resp = data - vmodp

            vmin = np.nanmin(data)
            vmax = np.nanmax(data)

            im = ax1.matshow(data, origin='lower', vmin=vmin, vmax=vmax)
            im = ax2.matshow(vmodp, origin='lower', vmin=vmin, vmax=vmax)
            im = ax3.matshow(resp, origin='lower', vmin=vmin, vmax=vmax)

            """
            im = ax3.matshow(resp, origin='lower',\
                        vmin=np.percentile(resp, 25),\
                vmax=np.percentile(resp, 75))
            if j == 0:
                cbar = fig.colorbar(im)
            else:
                cbar.draw_all()
            """

            plt.draw()

        rmin = rmaj
        rmaj = rmin + dr
        j += 1

    mask = mask.reshape(data.shape[0], data.shape[1])
    vmod = vmod.reshape(data.shape[0], data.shape[1])
    rotmod = rotmod.reshape(data.shape[0], data.shape[1])

    return rads, hds, ephds, emhds, mask, vmod, rotmod


if __name__ == "__main__":

    from galaxy_params import gal_feath as gal
    from constants import ang_to_phys

    # path = '/mnt/MyRAID/M33/VLA/14B-088/HI/full_imaging_noSD/'
    # filename = os.path.join(
    #     path, 'M33_14B-088_HI.clean.image.pbcov_gt_0.5_masked.mom1.fits')
    # out_path = os.path.join(path, "harmdec")
    # out_name = "M33_14B-088_HI.clean.image.pbcov_gt_0.5_masked.mom1"

    # path = '/mnt/MyRAID/M33/VLA/14B-088/HI/full_imaging_wGBT/'
    path = '/mnt/bigdata/ekoch/M33/VLA/14B-088/HI/full_imaging_wGBT/'
    filename = os.path.join(
        path, 'M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.mom1.fits')
    out_path = os.path.join(path, "harmdec")
    out_name = "M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.mom1"

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    beam = 18.9  # in arcsec
    dist = 840e3  # in pc
    # vsys = -179.56  # in km/s
    vsys = gal.vsys.to(u.km / u.s).value
    # pa = 200.66  # in deg
    pa = gal.position_angle.value  # in deg
    # incl = 58.  # in deg
    incl = gal.inclination.value  # in deg
    # ra0 = 23.4575 # of the center
    ra0 = gal.center_position.ra.value
    # dec0 = 30.6586  # of the center
    dec0 = gal.center_position.dec.value
    # ra0 = 15 * (1. + 33. / 60 + 50.9 / 3600)  # of the center
    # dec0 = 30. + 39. / 60 + 39. / 3600  # of the center
    rmax = 1980.  # in arcsec
    as2pc = (np.pi / 180) * np.asarray(dist) / 3600.

    mom1 = fits.open(filename)[0]
    data = mom1.data.squeeze()
    hd = mom1.header
    data[data == data[0, 0]] = np.nan

    data = data / 1000. - vsys  # everything in km/s with subtracted vsys

    order = 5  # harmonic decomposition order

    props = [ra0, dec0, dist, incl, pa, 0., rmax, beam]
    # rads, hds, ephds, emhds, mask, rec, rotmod = \
    #     harmdec(data, hd, props, order=order, blind=True)
    out = harmdec(data, hd, props, order=order, blind=True)

    rec[rec == 0] = np.nan
    print filename, "Residual = ", np.nanmedian(data - rec)

    # Harmonic terms table
    hdtab = Table([rads, hds[0], emhds[0], ephds[0]],
                  names=('Radius', 'c0', 'emc0', 'epc0'))

    for j in range(order):
        hdtab.add_column(Column(data=hds[j + 1], name='c' + str(j + 1)))
        hdtab.add_column(Column(data=emhds[j + 1], name='emc' + str(j + 1)))
        hdtab.add_column(Column(data=ephds[j + 1], name='epc' + str(j + 1)))
        hdtab.add_column(Column(data=hds[j + 2], name='s' + str(j + 1)))
        hdtab.add_column(Column(data=emhds[j + 2], name='ems' + str(j + 1)))
        hdtab.add_column(Column(data=ephds[j + 2], name='eps' + str(j + 1)))

    hdtab.write(os.path.join(out_path, out_name + '_hd.txt'), format='ascii',
                overwrite=True)

    # Save the model images
    mask = fits.PrimaryHDU(mask, hd)
    mask.writeto(os.path.join(out_path, out_name + '_mask.fits'),
                 overwrite=True)
    rec = fits.PrimaryHDU(rec, hd)
    rec.writeto(os.path.join(out_path, out_name + '_rec.fits'),
                overwrite=True)
    rotmod = fits.PrimaryHDU(rotmod, hd)
    rotmod.writeto(os.path.join(out_path, out_name + '_rotmod.fits'),
                   overwrite=True)
