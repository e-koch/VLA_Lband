import os
import sys
import time
import pickle

from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
from imagerhelpers.input_parameters import ImagerParameters

# CASA init should have the VLA_Lband repo appended to the path
from paths import data_path

full_path = os.path.join(data_path, "14B-088")

scratch_path = sys.argv[-1]
output_path = os.path.join(scratch_path, "dirty_cube")

if not os.path.exists(output_path):
    os.mkdir(output_path)

orig_dir = os.getcwd()

# Change out to the scratch directory
os.chdir(scratch_path)

paramList = \
    ImagerParameters(msname='14B-088_HI.ms.contsub',
                     datacolumn='data',
                     field='M33*',
                     imagename='dirty_cube/M33_14B-088_HI.dirty',
                     # imagename=os.path.join(output_path[-1],
                     #                        'M33_14B-088_HI.dirty'),
                     imsize=[2560, 2560],
                     cell='3arcsec',
                     specmode='cube',
                     start=800,
                     width=1,
                     nchan=64,
                     startmodel=None,
                     gridder='mosaic',
                     weighting='natural',
                     niter=0,
                     threshold='3.2mJy/beam',
                     phasecenter='J2000 01h33m50.904 +30d39m35.79',
                     restfreq='1420.40575177MHz',
                     outframe='LSRK',
                     pblimit=0.2,
                     usemask='pb',
                     mask=None,
                     deconvolver='hogbom',
                     dopbcorr=False,
                     chanchunks=-1
                     )

imager = PyParallelCubeSynthesisImager(params=paramList)

# init major cycle elements
imager.initializeImagers()
imager.initializeNormalizers()
imager.setWeighting()

# Init minor cycle elements
imager.initializeDeconvolvers()
imager.initializeIterationControl()
imager.makePSF()
imager.makePB()

# Make dirty image
t0 = time.time()
imager.runMajorCycle()
t1 = time.time()
casalog.post("Time for major cycle: {}".format(t1 - t0))

retrec = imager.getSummary()

concattype = 'virtualcopy'
imager.concatImages(type=concattype)
imager.deleteTools()

# Save the output dict
with open(os.path.join(output_path,
                       'M33_14B-088_HI.dirty.clean_output.pickle'),
          'wb') as handle:
    pickle.dump(retrec, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Return to original directory
os.chdir(orig_dir)
