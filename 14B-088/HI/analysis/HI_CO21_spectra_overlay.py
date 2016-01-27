
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.pyplot as p

'''
Overlay the CO (2-1) IRAM spectra w/ an HI spectrum with absorption features
'''


hi_cube = SpectralCube.read(
    '/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.fits',
    mode='denywrite')

hi_old_cube = SpectralCube.read(
    '/media/eric/Data_3/M33/VLA_Data/AT0206/imaging/M33_206_b_c_HI.fits',
    mode='denywrite')

# arec_cube = "/media/eric/Data_3/M33/Arecibo/14B-088_items/M33_14B-088_HI_model.fits"
# arec_cube = "/media/eric/Data_3/M33/Arecibo/M33only.fits"

# arecibo = SpectralCube.read(arec_cube, mode='denywrite')

co21_cube = SpectralCube.read('/media/eric/Data_3/M33/IRAM/m33.co21_iram.fits')
co21_cube = co21_cube.spectral_slab(*hi_cube.spectral_extrema)

# Only keep the overlap
lat_extrema = co21_cube.latitude_extrema
long_extrema = co21_cube.longitude_extrema
hi_cube = hi_cube.subcube(xlo=long_extrema[1],
                          xhi=long_extrema[0],
                          ylo=lat_extrema[0],
                          yhi=lat_extrema[1])

# Now find the position in each cube corresponding to a given
# coordinate

# posn = SkyCoord("01h33m21.287", "+30d32m16.110", frame='icrs')
# posn = SkyCoord("01h34m31.492", "+30d46m32.248", frame='icrs')
# posn = SkyCoord("01h34m31.87", "+30d46m56.5", frame='icrs')
posn = SkyCoord("01h34m34.65", "+30d46m33.73", frame='icrs')


def get_closest_posn(posn, spatial_footprint):
    '''
    '''

    spatial_coords = SkyCoord(ra=spatial_footprint[1],
                              dec=spatial_footprint[0],
                              frame='icrs')

    min_posn = spatial_coords.separation(posn).argmin()

    twod_posn = np.unravel_index(min_posn, spatial_footprint[0].shape)

    return twod_posn


def plot_overlay_spectra(spectrum1, spectrum2, figsize=(12, 4),
                         spec_unit=u.km/u.s, label1="HI", label2="CO (2-1)",
                         style='steps-mid', xlim=None, second_axis=True,
                         ylabel="HI Surface Brightness", vline=None,
                         hline=None):
    '''
    Overlay 2 spectra with individual y-axes.

    Parameters
    ----------
    vline : float or 3-element tuple.
        If a tuple is given, the 2nd and 3rd values for the lower and upper
        bounds.
    '''

    p.figure(figsize=figsize)
    ax = p.subplot(111)
    ax.plot(spectrum1.spectral_axis.to(spec_unit), spectrum1.value,
            drawstyle=style, label=label1)
    ax.set_xlabel('Velocity ('+str(spec_unit)+')')
    if second_axis:
        ax.set_ylabel(label1+' ('+str(spectrum1.unit)+')')
    else:
        ax.set_ylabel(ylabel+' ('+str(spectrum1.unit)+')')
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

    if second_axis:
        ax_2 = ax.twinx()
        ax_2.plot(spectrum2.spectral_axis.to(spec_unit), spectrum2.value, 'g',
                  drawstyle=style)
        ax_2.set_ylabel(label2+' ('+str(spectrum2.unit)+')')
        align_yaxis(ax, 0, ax_2, 0)
    else:
        ax.plot(spectrum2.spectral_axis.to(spec_unit), spectrum2.value, 'g',
                drawstyle=style, label=label2)
        ax.legend()

    if xlim is not None:
        ax.set_xlim(xlim)

    if vline is not None:
        ymin, ymax = ax.get_ylim()
        ax.vlines(vline[0], ymin, ymax, colors='r',
                  linestyles="--")

        diff = np.abs(spectrum1.spectral_axis[1] -
                      spectrum1.spectral_axis[0]).to(spec_unit).value

        print diff
        ax.fill_between(np.arange(vline[1], vline[2]+diff, diff),
                        ymin, ymax, facecolor='r', alpha=0.25)

    if hline is not None:
        ax.hlines(hline, ax.get_ylim()[0], ax.get_ylim()[1], colors='b',
                  linestyles="--")

    p.tight_layout()

# Helper functions from
# http://stackoverflow.com/questions/10481990/matplotlib-axis-with-two-scales-shared-origin


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2, (y1-y2)/2, v2)
    adjust_yaxis(ax1, (y2-y1)/2, v1)


def adjust_yaxis(ax, ydif, v):
    """shift axis ax by ydiff, maintaining point v at the same location"""
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny > maxy or (-miny == maxy and dy > 0):
        nminy = miny
        nmaxy = miny*(maxy+dy)/(miny+dy)
    else:
        nmaxy = maxy
        nminy = maxy*(miny+dy)/(maxy+dy)
    ax.set_ylim(nminy+v, nmaxy+v)

hi_posn = get_closest_posn(posn, hi_cube.spatial_coordinate_map)
hi_old_posn = get_closest_posn(posn, hi_old_cube.spatial_coordinate_map)

co21_posn = get_closest_posn(posn, co21_cube.spatial_coordinate_map)

# arec_posn = get_closest_posn(posn, arecibo.spatial_coordinate_map)

# Convert to K
hi_spectrum = hi_cube[:, hi_posn[0], hi_posn[1]]
hi_spectrum = hi_spectrum.to(u.K, hi_cube.beam.jtok_equiv(1420.40575177*u.MHz))

hi_old_spectrum = hi_old_cube[:, hi_old_posn[0], hi_old_posn[1]]
hi_old_spectrum = hi_old_spectrum.to(u.K, hi_old_cube.beam.jtok_equiv(1420.40575177*u.MHz))

# hi_arec_spectrum = arecibo[:, arec_posn[0], arec_posn[1]]
# hi_arec_spectrum = hi_arec_spectrum.to(u.K, arecibo.beam.jtok_equiv(1420.40575177*u.MHz))

co21_spectrum = co21_cube[:, co21_posn[0], co21_posn[1]]

# hi_spectrum.quicklook(label='HI', color='r')
# co21_spectrum.quicklook(label='CO (2-1)', color='b')

plot_overlay_spectra(hi_spectrum, hi_old_spectrum,
                     label1="EVLA",
                     label2="Historical VLA", second_axis=False,
                     vline=(-226, -238, -214))

plot_overlay_spectra(hi_spectrum, co21_spectrum,
                     label1="HI Brightness",
                     label2="CO(2-1) Brightness",
                     vline=(-226, -238, -214))
