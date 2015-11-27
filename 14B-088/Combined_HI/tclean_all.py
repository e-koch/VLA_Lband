
'''
Using tclean to image the regridded 14B-088 HI data and the archival data
together.
'''

from tasks import tclean

vis = ["14B-088_HI_LSRK_AT0206_regrid.ms.contsub",
       "M33_b_c_LSRK.ms"]
imagename = "M33_AT0206_14B-088_HI.clean"
model = "M33_model.image"
mask = "M33_newmask.image"

tclean(vis=vis, imagename=imagename, field="M33*", imsize=[4096, 4096],
       cell="1.5arcsec", phasecenter="J2000 01h33m50.904 +30d39m35.79",
       restfreq="1420.40575177MHz", startmodel=model, specmode='cube',
       nchan=255, start=1, width=1, veltype='radio', gridder='mosaic',
       pblimit=0.2, deconvolver="multiscale", scales=[0, 4, 8, 20, 40, 80],
       weighting='natural', niter=80000, threshold="2mJy/beam",
       cycleniter=5000, mask=mask, parallel=True)
