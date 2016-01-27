
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.wcs_utils import drop_axis
from pvextractor import Path, extract_pv_slice
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from aplpy import FITSFigure
import numpy as np
import matplotlib.pyplot as p

'''
Create a PV Slice
'''

# hi_cube = SpectralCube.read(
#     '/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.fits',
#     mode='denywrite')

# hi_cube = SpectralCube.read(
#     '/media/eric/Data_3/M33/VLA_Data/AT0206/imaging/M33_206_b_c_HI.fits',
#     mode='denywrite')

hi_cube = SpectralCube.read(
    '/media/eric/MyRAID/M33/14B-088/HI/combined_HI/M33_14B-088_AT0206_HI.clean.image.fits',
    mode='denywrite')

co10_cube = SpectralCube.read('/media/eric/Data_3/M33/co/m33.merge.20.fits')
co10_cube = co10_cube.spectral_slab(*hi_cube.spectral_extrema)

co21_cube = SpectralCube.read('/media/eric/Data_3/M33/co21/m33.co21_iram.fits')
co21_cube = co21_cube.spectral_slab(*hi_cube.spectral_extrema)

co32_cube = SpectralCube.read('/media/eric/Data_3/M33/co/m33.co32.fits')
co32_cube = co32_cube.spectral_slab(*hi_cube.spectral_extrema)

# Only keep the overlap
lat_extrema = co21_cube.latitude_extrema
long_extrema = co21_cube.longitude_extrema
hi_cube = hi_cube.subcube(xlo=long_extrema[1],
                          xhi=long_extrema[0],
                          ylo=lat_extrema[0],
                          yhi=lat_extrema[1])
hi_cube = hi_cube.spectral_slab(-210*u.km/u.s, -270*u.km/u.s)
co21_cube = co21_cube.spectral_slab(-210*u.km/u.s, -270*u.km/u.s)
co32_cube = co32_cube.spectral_slab(-210*u.km/u.s, -270*u.km/u.s)
co10_cube = co10_cube.spectral_slab(-210*u.km/u.s, -270*u.km/u.s)

# Simple masks
co21_sigma = np.nanstd(co21_cube[-1].value)
co32_sigma = np.nanstd(co32_cube[-1].value)
co10_sigma = np.nanstd(co10_cube[0].value)
# co21_cube = co21_cube.with_mask(co21_cube > co21_sigma*co21_cube.unit)


slice_centers = [SkyCoord("1h34m34.7", "+30d45m38.0"),
                 SkyCoord("1h34m34.7", "+30d45m55.0"),
                 SkyCoord("1h34m34.7", "+30d46m12.0"),
                 SkyCoord("1h34m34.7", "+30d46m29.0"),
                 SkyCoord("1h34m34.7", "+30d46m46.0"),
                 SkyCoord("1h34m34.7", "+30d47m03.0"),
                 SkyCoord("1h34m34.7", "+30d47m20.0"),
                 SkyCoord("1h34m34.7", "+30d47m37.0"),
                 SkyCoord("1h34m34.7", "+30d47m54.0"),
                 SkyCoord("1h34m34.7", "+30d48m11.0"),
                 ]

paths = [Path(SkyCoord([coord.ra + 80*u.arcsec, coord.ra - 80*u.arcsec],
                       coord.dec), width=17*u.arcsec)
         for coord in slice_centers]

# for path in paths:
#     pv = extract_pv_slice(hi_cube, path)
#     pv_co21 = extract_pv_slice(co21_cube, path)
#     pv_co32 = extract_pv_slice(co32_cube, path)
#     pv_co10 = extract_pv_slice(co10_cube, path)

#     fig = FITSFigure(pv)
#     fig.show_grayscale()
#     fig.show_colorbar()
#     fig.show_contour(pv_co21, levels=np.arange(3*co21_sigma, 10*co21_sigma,
#                                                2*co21_sigma),
#                      colors='g')
#     fig.show_contour(pv_co32, levels=np.arange(2*co32_sigma, 10*co32_sigma,
#                                                2*co32_sigma),
#                      colors='r')
#     fig.show_contour(pv_co10, levels=np.arange(2*co10_sigma, 10*co10_sigma,
#                                                2*co10_sigma),
#                      colors='b')
#     raw_input("?")
#     fig.close()

# Combine some of the interesting PV slices w/ moment0 maps:


def load_hi():
    hi_mom0 = \
        fits.open("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"
                  "M33_14B-088_HI.clean.image.mom0_3sigmask.fits")[0]

    hi_proj = Projection(hi_mom0.data, unit=u.Jy, wcs=WCS(hi_mom0.header))

    return hi_mom0, hi_proj


def load_halp_hubble():

    halpha = fits.open("/media/eric/Data_3/M33/Hubble/hst_05773_03_wfpc2_f656n_wf"
                       "/hst_05773_03_wfpc2_f656n_wf_sci.fits")[1]

    halpha_proj = Projection(halpha.data, unit=u.count,
                             wcs=WCS(halpha.header))

    return halpha, halpha_proj


def load_co21_IRAM():

    co21 = fits.open("/media/eric/Data_3/M33/co21/m33.ico.fits")[0]

    new_wcs = drop_axis(drop_axis(WCS(co21.header), 3), 2)
    co21_proj = Projection(co21.data.squeeze(), unit=u.K, wcs=new_wcs)

    return co21, co21_proj

# Cut down HI and CO21 to have the same coverage as Halpha
_, halp_proj = load_halp_hubble()
_, hi_proj = load_hi()
_, co21_proj = load_co21_IRAM()

halp_long = halp_proj.longitude_extrema
halp_lat = halp_proj.latitude_extrema

hi_proj = hi_proj.subprojection(xlo=halp_long[1], xhi=halp_long[0],
                                ylo=halp_lat[0], yhi=halp_lat[1])

co21_proj = co21_proj.subprojection(xlo=halp_long[1], xhi=halp_long[0],
                                    ylo=halp_lat[0], yhi=halp_lat[1])

fig = p.figure()

# Halpha
halp_fig = FITSFigure(halp_proj.hdu, subplot=(2, 3, 1), figure=fig)
halp_fig.show_grayscale()
halp_fig.hide_xaxis_label()
halp_fig.hide_yaxis_label()
halp_fig.hide_tick_labels()
halp_fig.set_title(r"H$\alpha$")

hi_fig = FITSFigure(hi_proj.hdu, subplot=(2, 3, 2), figure=fig)
hi_fig.show_grayscale()
hi_fig.hide_xaxis_label()
hi_fig.hide_yaxis_label()
hi_fig.hide_tick_labels()
hi_fig.set_title(r"HI")
hi_fig.set_nan_color("k")
hi_fig.show_regions('pvslices_finalfig.reg')

co21_fig = FITSFigure(co21_proj.hdu, subplot=(2, 3, 3), figure=fig)
co21_fig.show_grayscale()
co21_fig.hide_xaxis_label()
co21_fig.hide_yaxis_label()
co21_fig.hide_tick_labels()
co21_fig.set_title(r"CO(2-1)")

for path, pnum in zip(paths[4:7], [4, 5, 6]):
        pv = extract_pv_slice(hi_cube, path)
        pv_co21 = extract_pv_slice(co21_cube, path)

        pv_fig = FITSFigure(pv, subplot=(2, 3, pnum), figure=fig)
        pv_fig.show_grayscale()
        pv_fig.show_contour(pv_co21,
                            levels=np.arange(3*co21_sigma, 20*co21_sigma,
                                             4*co21_sigma),
                            colors='g')

        pv_fig.hide_xaxis_label()
        pv_fig.hide_xtick_labels()

        if pnum != 4:
            pv_fig.hide_yaxis_label()
            pv_fig.hide_ytick_labels()

p.tight_layout()
