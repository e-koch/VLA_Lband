
'''
Use SPW 3 and 6 from two tracks to make a nice 604 image.
Aim to combine with B and A configs over same region.
'''

myvis = ['14B-088_10_21_14/14B-088_10_21_14_continuum/14B-088.sb29701604.eb29882607.56952.08797296297.continuum.ms',
         '14B-088_11_10_14/14B-088_11_10_14_continuum/14B-088.sb29971045.eb29984617.56972.01955872685.continuum.ms']

# Setup flanking fields
# Cross-referenced bright features in NVSS and the full M33 mosaic
# 'J2000 1h32m22.618 +30d44m08.711'
# 'J2000 1h35m22.961 +30d27m29.49'
# 'J2000 1h37m26.8 +31d02m13.4'
# 'J2000 1h36m17.6 +31d24m32.7'
# 'J2000 1h36m54.2 +31d14m31.6'
# 'J2000 1h37m08.7 +31d22m36.3'

tclean(vis=myvis,
       datacolumn='corrected',
       imagename='M33_14B-088_NGC604_continuum_spw_3_6_big',
       field='M33_2,M33_3,M33_7_center,M33_8',
       spw="3:20~122,6",
       imsize=[5000, 5000],
       cell='3arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.000001,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48],
       pbcor=False,
       interactive=True,
       usemask='auto-multithresh',
       pbmask=0.0,
       # outlierfile='ngc604_outliers.txt',
       calcpsf=True,
       calcres=True,
       )


# Removed from outlier list

# imagename=outlier3
# imsize=[64, 64]
# phasecenter=J2000 01:37:26.8 +31.02.13.4
# deconvolver=hogbom
# gridder=standard

# imagename=outlier4
# imsize=[64, 64]
# phasecenter=J2000 01:36:17.6 +31.24.32.7
# deconvolver=hogbom
# gridder=standard

# imagename=outlier5
# imsize=[64,64]
# phasecenter=J2000 01:36:54.2 +31.14.31.6
# deconvolver=hogbom
# gridder=standard

# imagename=outlier6
# imsize=[64, 64]
# phasecenter=J2000 01:37:08.7 +31.22.36.3
# deconvolver=hogbom
# gridder=standard

# Try a outlier field approach with a single pointing

tclean(vis=myvis,
       datacolumn='corrected',
       imagename='M33_14B-088_NGC604_continuum_spw_3_6_ptg_8_only',
       field='M33_3',
       spw="3:20~122,6",
       imsize=[1000, 1000],
       cell='3arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='standard',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       # phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.000000001,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48],
       pbcor=False,
       interactive=True,
       usemask='auto-multithresh',
       pbmask=0.000000001,
       outlierfile='ngc604_outliers.txt',
       calcpsf=True,
       calcres=True,
       )


tclean(vis=myvis,
       datacolumn='corrected',
       imagename='outlier5_alone',
       field='M33_4',
       spw="3:20~122,6",
       imsize=[64, 64],
       cell='3arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='standard',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01:36:54.2 +31.14.31.6',
       # phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.0000001,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48],
       pbcor=False,
       interactive=True,
       usemask='pb',
       pbmask=0.0000001,
       calcpsf=True,
       calcres=True,
       )
