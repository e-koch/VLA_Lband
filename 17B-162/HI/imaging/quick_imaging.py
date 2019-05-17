
'''
Make a single pass, 17B and 14B cube, cleaned to 5-sigma.
Useful as a check against the old version while (hopefully)
taking less user-time to compared to the final version.
'''

import os

from tasks import tclean

myvis = ['17B-162_HI_spw_0_LSRK.ms.contsub',
         '14B-088_HI_LSRK.ms.contsub']

myfields = 'M33_2,M33_6,M33_7_center,M33_8,M33_11,M33_12,M33_14'

spw_num = 0

output_path = "HI_quick_imaging"

imgname = 'M33_HI_14B_17B_quickcube'

myimagesize = [6250, 6144]
mycellsize = '1.0arcsec'
mypblimit = 0.2

# Capture all emission near centre of galaxy
start_vel = "-110km/s"
width_vel = "1km/s"
nchans = 130

if not os.path.exists(output_path):
    os.mkdir(output_path)

default('tclean')

# Make a dirty cube

# if it already exists, skip to the cleaning call below

dirty_resid_name = os.path.join(output_path, "{}.residual".format(imgname))

if not os.path.exists(dirty_resid_name):

    casalog.post("Making dirty cube.")

    tclean(vis=myvis,
           datacolumn='corrected',
           imagename=os.path.join(output_path, imgname),
           spw=str(spw_num),
           field=myfields,
           imsize=myimagesize,
           cell=mycellsize,
           specmode='cube',
           start=start_vel,
           width=width_vel,
           nchan=nchans,
           startmodel=None,
           gridder='mosaic',
           weighting='natural',
           niter=0,
           threshold='',
           nsigma=5.,
           phasecenter="J2000 01h33m50.904 +30d39m35.79",
           restfreq="1.420405752GHz",
           outframe='LSRK',
           pblimit=mypblimit,
           deconvolver='multiscale',
           pbcor=False,
           veltype='radio',
           chanchunks=-1,
           restoration=False,
           parallel=True,
           usemask='auto-multithresh',
           mask=None,
           pbmask=0.1,
           verbose=True,
           )

# Copy the dirty residual into its own directory
dirty_cube_path = os.path.join(output_path, "dirty_cube")

if not os.path.exists(dirty_cube_path):
    os.mkdir(dirty_cube_path)

os.system("cp -r {0}.residual {1}".format(os.path.join(output_path, imgname),
                                          dirty_cube_path))

# Initially clean to 5-sigma WITH automasking!
casalog.post("Cleaning cube to 5-sigma.")

tclean(vis=myvis,
       datacolumn='corrected',
       imagename=os.path.join(output_path, imgname),
       spw=str(spw_num),
       field=myfields,
       imsize=myimagesize,
       cell=mycellsize,
       specmode='cube',
       start=start_vel,
       width=width_vel,
       nchan=nchans,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=1000000,
       threshold='',
       nsigma=5.,
       phasecenter="J2000 01h33m50.904 +30d39m35.79",
       restfreq="1.420405752GHz",
       outframe='LSRK',
       pblimit=mypblimit,
       deconvolver='multiscale',
       scales=[0, 6, 12, 24, 48, 96],
       pbcor=False,
       veltype='radio',
       chanchunks=-1,
       restoration=True,
       parallel=True,
       cycleniter=10000,  # Force a lot of major cycles
       usemask='auto-multithresh',
       mask=None,
       pbmask=0.1,
       minpercentchange=2.,
       noisethreshold=5.,
       lownoisethreshold=2.,
       sidelobethreshold=2.5,
       minbeamfrac=0.2,
       verbose=True,
       calcres=False,
       calcpsf=False,
       fastnoise=False,  # Use noise calc more robust for extended emission
       )
