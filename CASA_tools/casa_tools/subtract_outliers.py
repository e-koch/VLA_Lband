
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

from taskinit import tb
from tasks import split, uvsub, concat, clean, rmtables

from .graceful_error_catch import catch_fail


def subtract_outliers(vis, outlier_coords, field='M33*', split_fields=True,
                      stokes='I', interactive=True, weighting='natural',
                      threshold='5mJy/beam', cell='3arcsec', cleanup=False,
                      datacolumn="CORRECTED", imsize=64, save_space=False,
                      masks=None):
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

    if masks is not None:
        if len(masks) != len(outlier_coords):
            raise Warning("The number of specified masks must match the"
                          " number of coordinates.")

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

        # Image each outlier at its phasecenter, then uvsub
        for i, coord in enumerate(outlier_coords):

            outfield_img = fieldimg + "_" + str(i)

            if masks is not None:
                mask = masks[i]
            else:
                mask = None

            catch_fail(clean, vis=fieldvis, imagename=outfield_img, mode='mfs',
                       phasecenter=coord, niter=10000, usescratch=True,
                       interactive=interactive, cell='3arcsec',
                       imsize=imsize, threshold=threshold, weighting=weighting,
                       minpb=0.0, mask=mask)

            # Subtract out the model from the imaging
            catch_fail(uvsub, vis=fieldvis)

            # Remove the individual images
            if cleanup:
                catch_fail(rmtables, tablenames=outfield_img+"*")

        if save_space:

            fieldvis_corronly = \
                os.path.join('temp_files', f+"_corrected.ms")

            catch_fail(split, vis=fieldvis, outputvis=fieldvis_corronly,
                       field=f, datacolumn="CORRECTED")

            catch_fail(rmtables, tablenames=fieldvis)

    # Now append the uvsub fields back together
    individ_ms = glob.glob("temp_files/*.ms")
    catch_fail(concat, vis=individ_ms,
               concatvis=vis.rstrip(".ms")+"_outsub.ms",
               respectname=True)

    if cleanup:
        catch_fail(rmtables, tablenames='temp_files/*')
