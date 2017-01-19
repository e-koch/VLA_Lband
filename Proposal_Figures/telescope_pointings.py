
'''
Calculate a hexagonal grid of telescope pointings
'''

import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import astropy.units as u


def calc_pointings(obj_coord, orient=22, dist=0.277, grid=(4, 5, 4),
                   angle_offset=120 * u.deg):
    '''
    Given the central coords of an object, calculate a grid with a given
    spacing and grid.

    For a hexagonal grid, angle_offset should be 120 deg.

    Parameters
    ----------
    obj_coord : SkyCoord
        Center of the object.
    orient : float or astropy.Quantity
        Orientation of the semimajor axis relative to North.
    dist : float or astropy.Quantity
        Angle between grid points (beam radius). If a float is given, unit
        assumed to be degrees.
    grid : tuple
        Number of pointings in each row.

    Returns
    -------
    coord_grid : dict
        Dictionary of grid pointings. Keys correspond to the rows.
    '''

    # Each point is defined off of the object coordinates
    # If the number of points in the middle row is odd, the middle is the
    # object coords
    # Else, we calculate based on 1/2 distances of the given distance.

    try:
        assert orient.unit == u.deg
        assert dist.unit == u.deg
    except AttributeError:
        orient = orient * u.deg
        dist = dist * u.deg

    # Make dictionary for grid
    # Each key corresponds to a row
    coord_grid = dict()
    for i, row in enumerate(grid):
        coord_grid[i] = [None] * row

    if grid[1] % 2 != 0.0:
        mid = grid[1] / 2
        coord_grid[1][mid] = obj_coord

        for pt in range(1, (grid[1] / 2) + 1):
            coord_grid[1][mid - pt] = \
                sphere_point(obj_coord, orient, dist * pt)
            coord_grid[1][mid + pt] = \
                sphere_point(obj_coord,
                             orient + 180. * u.deg,
                             dist * pt)

    else:  # Define the two points closest to the center
        mid = grid[1] / 2 - 1  # Offset for 0th element

        for pt in range(1, (grid[1] / 2) + 1):
            coord_grid[1][mid - pt] = sphere_point(obj_coord, orient,
                                                   dist / 2. * pt)
            coord_grid[1][mid + pt] = sphere_point(obj_coord,
                                                   orient + 180 * u.deg,
                                                   dist / 2. * pt)

    # Now calculate the side rows based off of the centre row positions
    lower_center = sphere_point(obj_coord, orient + 90 * u.deg, dist)
    upper_center = sphere_point(obj_coord, orient - 90 * u.deg, dist)

    if grid[0] == 1:
        coord_grid[0][0] = lower_center
    else:
        for i in range(grid[0]):
            coord_grid[0][i] = sphere_point(coord_grid[1][i],
                                            orient + angle_offset,
                                            dist)
    if grid[2] == 1:
        coord_grid[2][0] = upper_center
    else:
        for i in range(grid[2]):
            coord_grid[2][i] = sphere_point(coord_grid[1][i],
                                            orient - angle_offset,
                                            dist)

    return coord_grid


def sphere_point(coord, orient, dist=0.277 * u.deg):
    '''
    Given a coordinate point and the bearing and distance, calculate
    the new point on the celestial sphere.

    Adapted the formula given at :
        http://trac.osgeo.org/openlayers/wiki/GreatCircleAlgorithms

    Parameters
    ----------
    coord : SkyCoord
        Reference coordinate.
    orient : float or astropy.Quantity
        Orientation of the direction of the new coordinate relative to North.
    dist : float or astropy.Quantity/ optional
        Angle between coordinates. If a float is given, unit
        assumed to be degrees.

    Returns
    -------
    new_coord : SkyCoord
        New Coordinates.
    '''

    # Check input units and tyoes
    assert isinstance(coord, SkyCoord)

    assert orient.unit == u.deg

    ra = coord.ra.to(u.deg)
    dec = coord.dec.to(u.deg)

    new_dec = np.arcsin(np.sin(dec.to(u.rad)) * np.cos(dist.to(u.rad)) +
                        np.cos(dec.to(u.rad)) * np.sin(dist.to(u.rad)) *
                        np.cos(orient.to(u.rad)))

    aa = np.sin(dist.to(u.rad)) * np.sin(orient.to(u.rad)) * u.rad

    bb = np.cos(dec.to(u.rad)) * np.cos(dist.to(u.rad)) - \
        np.sin(dec.to(u.rad)) * np.sin(dist.to(u.rad)) * np.cos(orient.to(u.rad))

    bb = bb * u.rad

    if bb != 0.0:
        new_ra = ra.to(u.deg) + np.arctan2(aa, bb).to(u.deg)
    else:
        new_ra = ra

    new_dec = new_dec.to(u.deg)

    new_coord = SkyCoord(new_ra, new_dec, 'icrs')

    assert np.isclose(coord.separation(new_coord).value, dist.value)

    return new_coord


if __name__ == "__main__":

    import regions as reg
    import os

    # Get the coords from Simbad
    m33 = Simbad.query_object("m33")
    ra = m33.field("RA").data[0]
    dec = m33.field("DEC").data[0]

    # Diameter of VLA dishes at 1.4 GHz
    dist = 0.277 * u.deg

    coord = SkyCoord(ra + " " + dec, 'icrs', unit=(u.hourangle, u.deg))

    # C config
    # grid = calc_pointings(coord, dist=dist.value)
    # B config
    grid = calc_pointings(coord, dist=dist.value, grid=(2, 3, 2),
                          angle_offset=120 * u.deg)
    # grid = calc_pointings(coord, dist=dist.value, grid=(1, 3, 1),
    #                       angle_offset=90 * u.deg)
    grid_3pt = calc_pointings(coord, dist=dist.value, grid=(0, 3, 0))

    for key in grid.keys():
        print "Row " + str(key)
        for posn in grid[key]:
            print(posn.to_string("hmsdms"))

    # Now convert into some regions
    regions = []
    for key in grid.keys():
        for posn in grid[key]:
            regions.append(reg.CircleSkyRegion(center=posn, radius=dist))

    regions_3pt = []
    for key in grid.keys():
        for posn in grid_3pt[key]:
            regions_3pt.append(reg.CircleSkyRegion(center=posn, radius=dist))

    # reg.writeds9(regions, "m33_cconfig.reg")
    out_path = \
        os.path.expanduser("~/Dropbox/code_development/VLA_Lband/Proposal_Figures")
    reg.write_ds9(regions, os.path.join(out_path, "m33_bconfig_7point.reg"))
    # reg.write_ds9(regions, "m33_bconfig_5point.reg")
    reg.write_ds9(regions_3pt, os.path.join(out_path, "m33_bconfig_3point.reg"))

