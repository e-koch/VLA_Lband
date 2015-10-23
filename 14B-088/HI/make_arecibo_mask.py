
'''
Create the mask with the Arecibo data, that is then interpolated to the VLA
cube and used as the CLEAN mask
'''

from astropy.io import fits


# I should really make this a package...
execfile("/home/eric/Dropbox/code_development/clean_masks/clean_mask_construct.py")

arecibo = \
    fits.open("/media/eric/Data_3/M33/Arecibo/M33only_jy_stokes_vrad.fits")[0]

cleanmask = CleanMask(arecibo, 1, 5, iteraxis=1)

# Need hard disable of slicewise std computation!!
cleanmask.make_mask(verbose=True, smooth=True, kern_size=6)

cleanmask.save_to_fits("M33only_jy_stokes_vrad_mask_1sigma_dilated.fits")
