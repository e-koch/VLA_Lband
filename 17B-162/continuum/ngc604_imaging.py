
'''
Use SPW 3 and 6 from two tracks to make a nice 604 image.
Aim to combine with C and A configs over same region.
'''

myvis = ['17B-162_09_23_17/17B-162_09_23_17_continuum/17B-162.sb34051874.eb34497862.58019.12101613426.continuum.ms',
         '17B-162_09_26_17/17B-162_09_26_17_continuum/17B-162.sb34051874.eb34512512.58022.09389369213.continuum.ms']

tclean(vis=myvis,
       datacolumn='corrected',
       imagename='M33_17B-162_NGC604_continuum_spw_3_6',
       field='M33_2,M33_7_center,M33_8',
       spw="3:20~122,6",
       imsize=[3000, 3000],
       cell='1arcsec',
       specmode='mfs',
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=10000,
       threshold='0.01mJy/beam',
       phasecenter='J2000 01h34m33.19 +30d47m05.6',
       pblimit=0.5,
       pbmask=0.5,
       deconvolver='hogbom',
       pbcor=False,
       interactive=True,
       usemask='auto-multithresh',
       )


# per-obs
for i, myvis_i in enumerate(myvis):

       tclean(vis=myvis_i,
              datacolumn='corrected',
              imagename='M33_17B-162_NGC604_continuum_spw_3_6_obs_{}'.format(i),
              field='M33_2,M33_7_center,M33_8',
              spw="3:20~122,6",
              imsize=[3000, 3000],
              cell='1arcsec',
              specmode='mfs',
              startmodel=None,
              gridder='mosaic',
              weighting='natural',
              niter=0,
              threshold='0.01mJy/beam',
              phasecenter='J2000 01h34m33.19 +30d47m05.6',
              pblimit=0.5,
              pbmask=0.5,
              deconvolver='hogbom',
              pbcor=False,
              interactive=False,
              usemask='auto-multithresh',
              restoration=False,
              )
