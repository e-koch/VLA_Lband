
'''
Create a catalogue of CO filament properties.
'''

from fil_finder import FilFinder2D
from spectral_cube import Projection
from astropy.io import fits
import astropy.units as u
import numpy as np
from scipy import ndimage as nd
from astropy.table import Table, Column

from paths import iram_co21_14B088_data_path, fourteenB_wGBT_HI_file_dict
from galaxy_params import gal_feath as gal

mom0 = Projection.from_hdu(fits.open(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.mom0.fits")))

hi_mom0 = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['Moment0']))

mask = nd.binary_erosion(np.isfinite(mom0), np.ones((3, 3)))
mask = nd.gaussian_filter(mask.astype(float), 1) >= 0.5

fil = FilFinder2D(mom0, distance=840 * u.kpc, mask=mask)
fil.preprocess_image(skip_flatten=False)
fil.create_mask(use_existing_mask=True,
                verbose=False)
fil.medskel(verbose=False)

fil.analyze_skeletons(prune_criteria='length',
                      branch_thresh=6 * u.pix,
                      skel_thresh=6 * u.pix)

# Create a table of skeleton properties

radii = gal.radius(header=hi_mom0.header).to(u.kpc)
pas = gal.position_angles(header=hi_mom0.header).to(u.deg)

dec, ra = mom0.spatial_coordinate_map

posn = np.where(fil.skeleton)

flat_ends = []
for end in fil.end_pts:
    for en in end:
        flat_ends.append(en)
flat_ends = np.array(flat_ends)

flat_inters = []
for inter in fil.intersec_pts:
    if len(inter) == 0:
        continue
    for ints in inter:
        for inn in ints:
            flat_inters.extend(ints)
flat_inters = np.array(flat_inters)

end_posn = []
for end in flat_ends:
    x_posn = posn[0] == end[0]
    y_posn = posn[1] == end[1]

    pos = np.argwhere(np.logical_and(x_posn, y_posn))[0]

    end_posn.append(pos[0])

is_end = np.zeros_like(posn[0], dtype=bool)
is_end[np.array(end_posn)] = True

assert is_end.sum() == flat_ends.shape[0]

inter_posn = []
for inter in flat_inters:
    x_posn = posn[0] == inter[0]
    y_posn = posn[1] == inter[1]

    pos = np.argwhere(np.logical_and(x_posn, y_posn))[0]

    inter_posn.append(pos[0])

is_inter = np.zeros_like(posn[0], dtype=bool)
is_inter[np.array(inter_posn)] = True

# There were some duplicates, so this won't quite be equal
# assert is_inter.sum() == flat_inters.shape[0]

tab = Table([Column(posn[0]), Column(posn[1]),
             Column(dec[posn]), Column(ra[posn]),
             Column(radii[posn]), Column(pas[posn]),
             Column(mom0[posn]), Column(is_inter),
             Column(is_end)],
            names=['ypix', 'xpix', 'dec', 'ra',
                   'radius', 'posangle', 'intensity',
                   'is_inter', 'is_end'])

tab.write(iram_co21_14B088_data_path("tables/co_filament_table.csv",
                                     no_check=True),
          format='ascii.ecsv')
