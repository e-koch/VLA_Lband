
'''
Given the location of outlier(s), check each field in a mosaic for coverage of
that location, image it, subtract out of the UV plane
'''

import sys
import os
import warnings
import re
import glob
from astropy.extern import six

from .graceful_error_catch import catch_fail


def subtract_outlier(vis, outlier_coords, field='M33*', split_fields=True,
                     stokes='I', interactive=True, weighting='natural',
                     threshold='5mJy/beam', cell='3arcsec', cleanup=False,
                     datacolumn="CORRECTED"):
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

    # If only one position given, convert to list for iteration.
    if isinstance(outlier_coords, six.string_types):
        outlier_coords = [outlier_coords]

    try:
        os.mkdir('temp_files')
    except OSError:
        warnings.warn("temp_files already exists. "
                      "Going to remove all ms files from it.")
        catch_fail(rmtables, tablenames='temp_files/*')

    # Loop through the fields
    for f in fields:

        fieldvis = os.path.join('temp_files', f+".ms")
        fieldimg = os.path.join('temp_files', f)

        # Split the field off
        catch_fail(split, vis=vis, outputvis=fieldvis, field=f,
                   datacolumn=datacolumn)

        # Now image the data, keeping the phasecenter at the outlier

        if interactive:
            mask = None
        else:
            # Make a mask at the center with a radius
            # slightly larger than the beam
            # NOT GENERALIZED!!
            mask = 'circle [ [ 32pix , 32pix] ,8pix ]'

        # Image each outlier at its phasecenter, then uvsub
        for i, coord in enumerate(outlier_coords):

            outfield_img = fieldimg + "_" + str(i)

            catch_fail(clean, vis=fieldvis, imagename=outfield_img, mode='mfs',
                       phasecenter=coord, niter=10000, usescratch=True,
                       interactive=interactive, cell='3arcsec',
                       imsize=64, threshold=threshold, weighting=weighting,
                       minpb=0.0)

            # Subtract out the model from the imaging
            catch_fail(uvsub, vis=fieldvis)

            # Remove the individual images
            if cleanup:
                catch_fail(rmtables, tablenames=outfield_img+"*")

    # Now append the uvsub fields back together
    individ_ms = glob.glob("temp_files/*.ms")
    catch_fail(concat, vis=individ_ms,
               concatvis=vis.rstrip(".ms")+"_outsub.ms",
               respectname=True)

    if cleanup:
        catch_fail(rmtables, tablenames='temp_files/*')
