
# Make the mask and the model

python ~/Dropbox/code_development/clean_masks/match_regrid_cubes.py /media/eric/Data_3/M33/Arecibo/M33only_jy_stokes_vrad.fits ~/MyRAID/M33/14B-088/HI/full_imaging/14B-088_HI.clean.image.fits /media/eric/Data_3/M33/Arecibo/M33_14B-088_HI_model.fits False

python ~/Dropbox/code_development/clean_masks/match_regrid_cubes.py /media/eric/Data_3/M33/Arecibo/M33only_jy_stokes_vrad_mask.fits ~/MyRAID/M33/14B-088/HI/full_imaging/14B-088_HI.clean.image.fits /media/eric/Data_3/M33/Arecibo/M33_14B-088_HI_mask.fits True

# Now swap the model and mask axes in the regridded cubes!