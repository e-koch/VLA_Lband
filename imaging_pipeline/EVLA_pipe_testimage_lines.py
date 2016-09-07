
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
if len(imaging_sources) == 0:
    print("No sources given. Input which field(s) should be imaged (mosaics"
          " can be created by giving a common name for the set; i.e., 'M33'"
          " for 'M33LP1', 'M33LP2', etc)")
    print("Multiple images can be created by separating the list w/ commas"
          " (i.e., '3C48, M33')")
    imaging_sources = raw_input("Input fields to be imaged: ")

    # Remove whitespaces then split by commas
    imaging_sources = imaging_sources.replace(" ", "").split(",")

# Make directory for images to go.
if not os.path.exists('test_images'):
    os.mkdir("test_images")

# Check list of given sources against the field list
valid_sources = []
for source in imaging_sources:
    if has_field(ms_active, source):
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

# Most common setups are going to have the same SPW coverage for each field
# So set spws to be the first one from field_spws
spws = field_spws[0]

for source in imaging_sources:
    for idx, spw_num in enumerate(spws):
        logprint("Imaging SPW {0} of {1}".format(idx, len(spws)),
                 logfileout='logs/testimage_lines.log')

        default("clean")

        weighting = 'natural'
        minpb = 0.1

        # Determine imagermode, cell size, and image size
        imagermode = set_imagermode(ms_active, source)
        cellsize = set_cellsize(ms_active, spw_num, sample_factor=6.)
        imagesize = set_imagesize(ms_active, spw_num, source, sample_factor=6.,
                                  pblevel=minpb)

        if imagermode == "mosaic":
            # XXX Set this to centre of M33 for now.
            phasecenter = 'J2000 01h33m50.904 +30d39m35.79'
        else:
            phasecenter = ''

        # Instead of trying to image a cube, just use one large channel
        # at the center (where large = 10 channels).
        center_chan = channels[idx] // 2
        width = 10
        start_chan = int(center_chan - width // 2)

        clean(vis=ms_active,
              imagename='test_images/{0}.{1}.spw_{2}'.\
                        format(ms_active[:-3], source, spw_num),
              field='*' + source + '*', spw=str(spw_num),
              mode='channel', niter=0,
              imagermode=imagermode, cell=cellsize,
              imsize=imagesize,
              start=start_chan, width=width, nchan=1,
              weighting=weighting, pbcor=False, minpb=0.1,
              phasecenter=phasecenter)


logprint("Finished EVLA_pipe_testimage_lines.py",
         logfileout='logs/testimage_lines.log')

pipeline_save()
