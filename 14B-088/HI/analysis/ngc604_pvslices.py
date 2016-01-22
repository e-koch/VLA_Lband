
from spectral_cube import SpectralCube
from pvextractor import Path, extract_pv_slice
from astropy.coordinates import SkyCoord
import astropy.units as u
from aplpy import FITSFigure
import numpy as np

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

for path in paths:
    pv = extract_pv_slice(hi_cube, path)
    pv_co21 = extract_pv_slice(co21_cube, path)
    pv_co32 = extract_pv_slice(co32_cube, path)
    pv_co10 = extract_pv_slice(co10_cube, path)

    fig = FITSFigure(pv)
    fig.show_grayscale()
    fig.show_colorbar()
    fig.show_contour(pv_co21, levels=np.arange(3*co21_sigma, 10*co21_sigma,
                                               2*co21_sigma),
                     colors='g')
    fig.show_contour(pv_co32, levels=np.arange(2*co32_sigma, 10*co32_sigma,
                                               2*co32_sigma),
                     colors='r')
    fig.show_contour(pv_co10, levels=np.arange(2*co10_sigma, 10*co10_sigma,
                                               2*co10_sigma),
                     colors='b')
    raw_input("?")
    fig.close()
