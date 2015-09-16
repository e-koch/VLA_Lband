
'''
Given the location of outlier(s), check each field in a mosaic for coverage of
that location, image it, subtract out of the UV plane
'''

import sys
import os
import warnings
import re
import glob


def subtract_outlier(vis, outlier_coord, field='M33*', split_fields=True,
                     stokes='I', interactive=True, weighting='natural',
                     threshold='5mJy/beam', cleanup=False):
    '''
    Subtract an outlier at the given coordinates. Splits out each field,
    tries to image at that coordinate, then subracts of the model in the UV
    plane. The fields are concatenated into a final ms.
    '''

    vis = vis.rstrip("/")

    # Get some info from the image
    tb.open(os.path.join(vis, 'FIELD'))
    all_fields = tb.getcol('NAME')
    tb.close()

    regex = re.compile(field)
    fields = [f for f in all_fields if re.match(regex, f)]

    try:
        os.mkdir('temp_files')
    except OSError:
        warnings.warn("temp_files already exists. "
                      "Going to remove all ms files from it.")
        rmtables('temp_files/*')

    # Loop through the fields
    for f in fields:

        fieldvis = os.path.join('temp_files', f+".ms")
        fieldimg = os.path.join('temp_files', f)

        # Split the field off
        split(vis=vis, outputvis=fieldvis, field=f,
              datacolumn="DATA")

        # Now image the data, keeping the phasecenter at the outlier

        if interactive:
            mask = None
        else:
            # Make a mask at the center with a radius
            # slightly larger than the beam
            # NOT GENERALIZED!!
            mask = 'circle [ [ 32pix , 32pix] ,8pix ]'

        clean(vis=fieldvis, imagename=fieldimg, mode='mfs',
              phasecenter=outlier_coord, niter=10000, usescratch=True,
              interactive=interactive, cell='3arcsec',
              imsize=64, threshold=threshold, weighting=weighting,
              minpb=0.0)

        # Subtract out the model from the imaging
        uvsub(vis=fieldvis)

    # Now append the uvsub fields back together
    individ_ms = glob.glob("temp_files/*.ms")
    concat(vis=individ_ms, concatvis=vis.rstrip(".ms")+"_outsub.ms",
           respectname=True)

    if cleanup:
        rmtables('temp_files/*')