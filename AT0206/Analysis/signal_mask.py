
'''
Make a mask of the emission
'''

cube = SpectralCube.read("M33_206_b_c_HI.fits")
cube = cube.with_mask(cube != 0*u.Jy)

noise_cube = Noise(cube)

new_noise = noise_cube.scale

cube = cube.with_mask(cube > new_noise*u.Jy)

# Load in the broad clean mask used
clean_mask = fits.getdata("../../../Arecibo/M33_newmask.fits")

# Need to match the dims
clean_mask = clean_mask.squeeze()
clean_mask = clean_mask[11:195, ::-1, ::-1]
clean_mask = clean_mask[:, 595:3504, 1065:3033]

old_mask = RadioMask(old_cube)
old_mask.intersection(old_arecibo_mask)
