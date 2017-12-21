
'''
Do large-scale OH stacking based on the peak HI velocity
'''

from spectral_cube.analysis_utilities import stack_spectra
from spectral_cube import SpectralCube, Projection
from astropy.io import fits
import astropy.units as u
import numpy as np

from paths import fourteenB_wGBT_HI_file_dict
from galaxy_params import gal_feath as gal

peakvels = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['PeakVels'])[0])

data_path = "/home/ekoch/bigdata/ekoch/VLA/14B-088/Lines/OH/"

oh_lines = ['OH1612', 'OH1665', 'OH1667', 'OH1720']

weighting = 'natural'

min_pb = 0.7

for line in oh_lines:

    cube_name = "{0}/{1}/imaging_1point5km_s/{1}_14B-088_{2}.image.pbcor.fits"\
        .format(data_path, line, weighting)
    pbcov_name = "{0}/{1}/imaging_1point5km_s/{1}_14B-088_{2}.pb.fits"\
        .format(data_path, line, weighting)

    oh_cube = SpectralCube.read(cube_name)
    oh_pb = fits.open(pbcov_name)[0]
    oh_cube = oh_cube.with_mask(oh_pb.data > min_pb)

    # Regrid the HI peak velocity to the OH grid
    # Already at a lower spatial resolution, so don't deconvolve
    reproj_hdr = oh_cube[0].wcs.to_header()
    reproj_hdr['NAXIS'] = 2
    reproj_hdr['NAXIS2'] = oh_cube.shape[1]
    reproj_hdr['NAXIS1'] = oh_cube.shape[2]
    peakvels_reproj = peakvels.reproject(reproj_hdr)

    radii = gal.radius(header=reproj_hdr)

    mask = radii <= 4 * u.kpc
    xy_posns = np.where(mask)

    stack_spec = stack_spectra(oh_cube, peakvels_reproj.quantity,
                               stack_function=np.nanmedian,
                               num_cores=1,
                               progressbar=True,
                               chunk_size=100000,
                               xy_posns=xy_posns)

    out_file = "{0}/{1}/imaging_1point5km_s/{1}_14B-088_{2}.HI_peak_stack.fits"\
        .format(data_path, line, weighting)

    stack_spec.write(out_file, overwrite=True)

    del oh_cube
