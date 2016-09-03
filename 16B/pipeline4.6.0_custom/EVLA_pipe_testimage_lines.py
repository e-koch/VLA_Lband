
'''
Create per SPW test images of a pipeline-calibrated MS.

This version is intended for line SPWs.

Run EVLA_pipe_restore.py before this script.
'''

logprint("Starting EVLA_pipe_testimage_lines.py",
         logfileout='logs/testimage_lines.log')

import os
import sys
from warnings import warn

from CASA_functions import (set_imagermode, set_imagesize, set_cellsize,
                            has_field)

# This script should still be usable if the user didn't enable imaging at the
# beginning. In this case, sources will be empty. Prompt the user at this
# point.
if len(sources) == 0:
    print("No sources given. Input which field(s) should be imaged (mosaics"
          " can be created by giving a common name for the set; i.e., 'M33'"
          " for 'M33LP1', 'M33LP2', etc)")
    print("Multiple images can be created by separating the list w/ commas"
          " (i.e., '3C48, M33')")
    sources = raw_input("Input fields to be imaged: ")
    # Remove whitespaces then split by commas
    sources = sources.replace(" ", "").split(",")

# Make directory for images to go.
if not os.path.exists('test_images'):
    os.mkdir("test_images")

# Check list of given sources against the field list
valid_sources = []
for source in sources:
    if has_field(vis, source):
        valid_sources.append(source)
    else:
        warn('No field contains the given source: {}'.format(source))
        logprint('No field contains the given source: {}'.format(source),
                 logfileout='logs/testimage_lines.log')

if len(valid_sources) == 0:
    warn('No valid sources given. Exiting without imaging.')
    logprint('No valid sources given. Exiting without imaging.',
             logfileout='logs/testimage_lines.log')
    sys.exit()

sources = valid_sources

for source in sources:
    for idx, spw_num in enumerate(spws):
        logprint("Imaging SPW {0} of {1}".format(idx, len(spws)),
                 logfileout='logs/testimage_lines.log')

        print("Imaging SPW {0} of {1}".format(idx, len(spws)))

        default("clean")
        # Determine imagermode
        imagermode = set_imagermode(vis, source)
        cellsize = set_cellsize(vis, spw_num, sample_factor=6.)
        imagesize = set_imagesize(vis, spw_num, sample_factor=6.)

        weighting = 'natural'
        # XXX Set this to centre of M33 for now.
        phasecenter = 'J2000 01h33m50.904 +30d39m35.79'

        # Instead of trying to image a cube, just use one large channel
        # at the center (where large = 10 channels).
        center_chan = channels[idx] // 2
        width = 10
        start_chan = center_chan - width // 2

        clean(vis=vis,
              imagename='test_images/{0}.{1}.spw_{2}'.format(vis[:-3], source,
                                                             spw_num),
              field=source + '*', spw=str(spw_num), mode='channel', niter=0,
              imagermode=imagermode, cell=cellsize, imagesize=imagesize,
              start=start_chan, width=width, nchan=1,
              weighting=weighting, pbcor=False, minpb=0.1,
              phasecenter=phasecenter)


logprint("Finished EVLA_pipe_testimage_lines.py",
         logfileout='logs/testimage_lines.log')

pipeline_save()
