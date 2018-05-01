
# Copy over the 14B-088 data to be released with the first M33 HI paper

getcert

DR1_PATH="vos:ekoch/M33/VLA/14B-088/DR1/"

DATA_PATH="/mnt/bigdata/ekoch/M33/VLA/14B-088/HI/full_imaging_wGBT"
DATA_PATH_noSD="/mnt/bigdata/ekoch/M33/VLA/14B-088/HI/full_imaging_noSD"

# HI Data products
DR1_HI_PATH=$DR1_PATH/HI/
vmkdir DR1_HI_PATH

# 2.6 km/s channel version of the cube
vcp $DATA_PATH/downsamp_to_co/M33_14B-088_HI.clean.image.GBT_feathered.2.6kms.fits $DR1_HI_PATH

# Primary beam coverage
vcp $DATA_PATH_noSD/M33_14B-088_pbcov.fits $DR1_HI_PATH

# Moment maps
vcp $DATA_PATH/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.kurtosis.fits $DR1_HI_PATH
vcp $DATA_PATH/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.lwidth.fits $DR1_HI_PATH
vcp $DATA_PATH/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.mom0.fits $DR1_HI_PATH
vcp $DATA_PATH/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.mom1.fits $DR1_HI_PATH
vcp $DATA_PATH/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peaktemps.fits $DR1_HI_PATH
vcp $DATA_PATH/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.peakvels.fits $DR1_HI_PATH
vcp $DATA_PATH/M33_14B-088_HI.clean.image.GBT_feathered.pbcov_gt_0.5_masked.skewness.fits $DR1_HI_PATH

# Stacked spectra
vmkdir $DR1_HI_PATH/stacked_spectra

vcp $DATA_PATH/stacked_spectra/centroid_stacked.fits $DR1_HI_PATH/stacked_spectra/
vcp $DATA_PATH/stacked_spectra/rotation_stacked.fits $DR1_HI_PATH/stacked_spectra/
vcp $DATA_PATH/stacked_spectra/peakvel_stacked.fits $DR1_HI_PATH/stacked_spectra/
vcp $DATA_PATH/stacked_spectra/*100pc*.fits $DR1_HI_PATH/stacked_spectra/
vcp $DATA_PATH/stacked_spectra/*5percentile*.fits $DR1_HI_PATH/stacked_spectra/
vcp $DATA_PATH/stacked_spectra/peak_stacking_5percentile_num_pix.npy $DR1_HI_PATH/stacked_spectra/

# Fits to the stacked spectra

vcp $DATA_PATH/tables/hi_hwhm_totalprof_fits_peak_5percentile_feather.csv $DR1_HI_PATH/stacked_spectra/

vcp $DATA_PATH_noSD/tables/hi_gaussian_hwhm_totalprof_fits_radial.csv $DR1_HI_PATH/stacked_spectra/
vcp $DATA_PATH_noSD/tables/hi_gaussian_totalprof_hwhm_fits.csv $DR1_HI_PATH/stacked_spectra/

# Rotation curve

vmkdir $DR1_HI_PATH/rotation_curve

vcp $DATA_PATH/diskfit_peakvels_noasymm_noradial_nowarp_output/diskfit_params_peakvels_nowarp_noradial_noasymm.inp $DR1_HI_PATH/rotation_curve/
vcp $DATA_PATH/diskfit_peakvels_noasymm_noradial_nowarp_output/rad.out* $DR1_HI_PATH/rotation_curve/
vcp $DATA_PATH/diskfit_peakvels_noasymm_noradial_nowarp_output/rad.fit*.*fits $DR1_HI_PATH/rotation_curve/
vcp $DATA_PATH/diskfit_peakvels_noasymm_noradial_nowarp_output/rad.mod.fits $DR1_HI_PATH/rotation_curve/
vcp $DATA_PATH/diskfit_peakvels_noasymm_noradial_nowarp_output/rad.res.fits $DR1_HI_PATH/rotation_curve/

# OH cubes

DR1_OH_PATH=$DR1_PATH/OH/
vmkdir $DR1_OH_PATH

OH_DATA_PATH="/mnt/bigdata/ekoch/VLA/14B-088/Lines/OH"

vcp $OH_DATA_PATH/OH1612/imaging_1point5km_s/OH1612_14B-088_natural.image.pbcor.fits $DR1_OH_PATH
vcp $OH_DATA_PATH/OH1665/imaging_1point5km_s/OH1665_14B-088_natural.image.pbcor.fits $DR1_OH_PATH
vcp $OH_DATA_PATH/OH1667/imaging_1point5km_s/OH1667_14B-088_natural.image.pbcor.fits $DR1_OH_PATH
vcp $OH_DATA_PATH/OH1720/imaging_1point5km_s/OH1720_14B-088_natural.image.pbcor.fits $DR1_OH_PATH


# RRL cubes

DR1_RRL_PATH=$DR1_PATH/RRL/
vmkdir $DR1_RRL_PATH

RRL_DATA_PATH="/mnt/bigdata/ekoch/VLA/14B-088/Lines/RRLs"

vcp $RRL_DATA_PATH/combined_imaging/HRL_14B-088_combined.fits $DR1_RRL_PATH

vcp $RRL_DATA_PATH/H152alp/dirty_images/H152alp_14B-088_dirty.image.fits $DR1_RRL_PATH
vcp $RRL_DATA_PATH/H152alp/dirty_images/H152alp_14B-088_dirty.image.smooth.fits $DR1_RRL_PATH

vcp $RRL_DATA_PATH/H153alp/dirty_images/H153alp_14B-088_dirty.image.fits $DR1_RRL_PATH
vcp $RRL_DATA_PATH/H153alp/dirty_images/H153alp_14B-088_dirty.image.smooth.fits $DR1_RRL_PATH

vcp $RRL_DATA_PATH/H158alp/dirty_images/H158alp_14B-088_dirty.image.fits $DR1_RRL_PATH
vcp $RRL_DATA_PATH/H158alp/dirty_images/H158alp_14B-088_dirty.image.smooth.fits $DR1_RRL_PATH

vcp $RRL_DATA_PATH/H164alp/dirty_images/H164alp_14B-088_dirty.image.fits $DR1_RRL_PATH
vcp $RRL_DATA_PATH/H164alp/dirty_images/H164alp_14B-088_dirty.image.smooth.fits $DR1_RRL_PATH

vcp $RRL_DATA_PATH/H166alp/dirty_images/H166alp_14B-088_dirty.image.fits $DR1_RRL_PATH
vcp $RRL_DATA_PATH/H166alp/dirty_images/H166alp_14B-088_dirty.image.smooth.fits $DR1_RRL_PATH

# Lastly copy over the README

vcp 14B-088/DR1_scripts/DR1_README.md $DR1_PATH
