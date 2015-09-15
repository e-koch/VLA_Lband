
'''
Make a mask of the emission
'''

from astropy.io import fits
from spectral_cube import SpectralCube, BooleanArrayMask
from signal_id import RadioMask, Noise
from astropy import units as u


make_mask = True
save_mask = False

cube = SpectralCube.read("M33_206_b_c_HI.fits")
cube = cube.with_mask(cube != 0*u.Jy)

if make_mask:

    # noise_cube = Noise(cube)

    # new_noise = noise_cube.scale

    new_noise = 0.0010592446454102172

    cube = cube.with_mask(cube > new_noise*u.Jy)

    # Load in the broad clean mask used
    clean_mask = fits.getdata("../../../Arecibo/M33_newmask.fits")

    # Need to match the dims
    clean_mask = clean_mask.squeeze()
    clean_mask = clean_mask[11:195, ::-1, ::-1]
    clean_mask = clean_mask[:, 595:3504, 1065:3033]

    from signal_id.utils import get_pixel_scales

    pixscale = get_pixel_scales(cube.wcs)

    beam_struct = cube.beam.as_tophat_kernel(pixscale)

    mask = RadioMask(cube)
    mask.intersection(clean_mask)
    # mask.remove_small_regions()
    # mask.open(iterations=3)
    # mask.close(iterations=3)
    # mask.dilate(iterations=6)


    if save_mask:
        mask.write('M33_206_b_c.source_mask.fits')
else:
    # Try loading the mask in

    mask_fits = fits.getdata('m33_206_b_c.source_mask.fits')

    mask = BooleanArrayMask(mask_fits.astype(bool), cube.wcs)

    cube = cube.with_mask(mask)
