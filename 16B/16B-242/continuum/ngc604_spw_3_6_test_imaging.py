
'''
Use first 4 tracks to image NGC 604 and the surrounding region.
'''

# myvis = ['16B-242_10_12_16/16B-242_10_12_16_continuum/16B-242.sb32681213.eb32949259.57673.06817918981.continuum.ms',
#          '16B-242_10_13_16/16B-242_10_13_16_continuum/16B-242.sb32681213.eb32954876.57674.07753513889.continuum.ms',
#          '16B-242_10_14_16/16B-242_10_14_16_continuum/16B-242.sb32681213.eb32956712.57675.0472325463.continuum.ms',
#          '16B-242_10_16_16/16B-242_10_16_16_continuum/16B-242.sb32681213.eb32958324.57677.085873668984.continuum.ms']

# tclean(vis=myvis,
#        datacolumn='corrected',
#        imagename='NGC604_16B-242_continuum_spw_3_6',
#        field='NGC604',
#        spw="3:20~122,6",
#        imsize=6000,
#        cell='0.25arcsec',
#        specmode='mfs',
#        startmodel=None,
#        gridder='mosaic',
#        weighting='natural',
#        niter=10000,
#        threshold='0.01mJy/beam',
#        phasecenter='J2000 01h34m33.19 +30d47m05.6',
#        pblimit=0.1,
#        pbmask=0.1,
#        deconvolver='hogbom',
#        pbcor=False,
#        interactive=True,
#        usemask='auto-multithresh',
#        )


# All configs

myvis = ['16B-242.sb32681213.eb32949259.57673.06817918981.continuum.ms',
         '16B-242.sb32681213.eb32954876.57674.07753513889.continuum.ms',
         '16B-242.sb32681213.eb32956712.57675.0472325463.continuum.ms',
         '16B-242.sb32681213.eb32958324.57677.085873668984.continuum.ms',
         '17B-162.sb34051874.eb34497862.58019.12101613426.continuum.ms',
         '17B-162.sb34051874.eb34512512.58022.09389369213.continuum.ms',
         '14B-088.sb29701604.eb29882607.56952.08797296297.continuum.ms',
         '14B-088.sb29971045.eb29984617.56972.01955872685.continuum.ms']

tclean(vis=myvis,
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_3_6_A_B_C',
       field='*',
       spw="3:20~122,6",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.01mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=True,
       usemask='auto-multithresh',
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=True,
       )


# SPW 3 only
tclean(vis=myvis,
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_3_A_B_C',
       field='*',
       spw="3:20~122",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=False,
       usemask='user',
       mask='NGC604_16B-242_continuum_spw_3_6_A_B_C.mask',
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=True,
       )

# SPW 6

tclean(vis=myvis,
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_6_A_B_C',
       field='*',
       spw="6",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=False,
       usemask='user',
       mask='NGC604_16B-242_continuum_spw_3_6_A_B_C.mask',
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=False,
       )

# A only

# Both SPWs with nterms=2
tclean(vis=myvis[:4],
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_3_6_A_mtmfs',
       field='*',
       spw="3:20~122,6",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       # deconvolver='multiscale',
       deconvolver='mtmfs',
       nterms=2,
       reffreq='1.5GHz',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=False,
       usemask='auto-multithresh',
       mask=None,
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=False,
       )


# SPW 3 only
tclean(vis=myvis[:4],
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_3_A',
       field='*',
       spw="3:20~122",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=False,
       usemask='user',
       mask='NGC604_16B-242_continuum_spw_3_6_A_B_C.mask',
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=False,
       )

# SPW 6

tclean(vis=myvis[:4],
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_6_A',
       field='*',
       spw="6",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=False,
       usemask='user',
       mask='NGC604_16B-242_continuum_spw_3_6_A_B_C.mask',
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=False,
       )


# All 4 still didn't give great dynamic range. Just try the first track?

tclean(vis=myvis[0],
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_3_A_1016',
       field='*',
       spw="3:20~122",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=False,
       usemask='auto-multithresh',
       mask=None,
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=False,
       )

# SPW 6

tclean(vis=myvis[0],
       datacolumn='corrected',
       imagename='NGC604_16B-242_continuum_spw_6_A_1016',
       field='*',
       spw="6",
       imsize=8000,
       cell='0.25arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.1mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.1,
       pbmask=0.1,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       interactive=False,
       usemask='auto-multithresh',
       mask=None,
       sidelobethreshold=3.0,
       noisethreshold=3.5,
       lownoisethreshold=1.5,
       negativethreshold=0.0,
       smoothfactor=2.0,
       minbeamfrac=0.1,
       cutthreshold=0.01,
       growiterations=75,
       fastnoise=False,
       parallel=False,
       )
