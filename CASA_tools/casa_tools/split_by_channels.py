
'''
Split MS into single spectral channels
'''

import os
import numpy as np
from copy import copy

from taskinit import tb, ia, rg

from .mytools import mymstransform, mysplit


def ms_split_by_channel(vis, nchan=-1, start=1, spw='0',
                        restfreq='1420.40575177MHz',
                        output_dir=None, use_split=True, **kwargs):
    '''
    Splits a MS by its spectral channels, according to the given
    '''

    # Get the total number of channels
    tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
    nchan_list = tb.getcol('NUM_CHAN')
    tb.close()

    # Pick out the spw specified
    if spw == "*":
        total_nchan = min(nchan_list)  # This isn't general at all...
    else:
        spw_slicer = get_slice_obj(spw)
        total_nchan = nchan_list[spw_slicer]

    # Now create a list of the channels to use
    if nchan == -1:
        nchan = total_nchan
    elif nchan > total_nchan:
        raise ValueError("There are only "+str(total_nchan)+". Reduce"
                         " nchan to at least this")
    else:
        pass

    channel_list = range(start, start+nchan)

    vis_name = vis.rstrip("/").split("/")[-1]

    for chan in channel_list:
        channel_vis = vis_name.rstrip(".ms")+"_channel_"+str(chan)+".ms"
        if output_dir is not None:
            channel_vis = os.path.join(output_dir, channel_vis)
        if use_split:
            spw_chan_select = spw + ":" + str(chan)
            mysplit(vis=vis, outputvis=channel_vis, spw=spw_chan_select,
                    **kwargs)
        else:
            mymstransform(vis=vis, outputvis=channel_vis, spw=spw,
                          regridms=True, width=1, nchan=nchan, start=chan,
                          restfreq=restfreq, **kwargs)


def image_split_by_channel(imagename, nchan=-1, start=1, output_dir=None,
                           specaxis_name="Frequency", verbose=False):
    '''
    Split an image into the specified spectral channels.
    '''

    # Open image.
    ia.open(imagename)

    # Find the spectral axis
    csys = ia.coordsys()
    try:
        spec_axis = np.where(np.asarray(csys.names()) == specaxis_name)[0][0]
    except IndexError:
        raise IndexError("Cannot find spectral axis" + specaxis_name + " in " +
                         str(csys.names()))
        ia.close()

    # Check given number of channels
    cube_shape = list(ia.shape())
    ndims = len(cube_shape)
    total_nchan = cube_shape[spec_axis]

    if nchan == -1:
        nchan = total_nchan
    elif nchan > total_nchan:
        raise ValueError("There are only "+str(total_nchan)+". Reduce"
                         " nchan to at least this")
        ia.close()
    else:
        pass

    # Create the base imagename
    output_dir = "" if output_dir is None else output_dir
    base_image = imagename.rstrip("/").split("/")[-1].rstrip(".image")

    outputimage = os.path.join(output_dir, base_image)

    channel_list = range(start, start+nchan)

    for chan in channel_list:
        if verbose:
            print("On channel "+str(chan+1)+" of "+str(start+nchan))
        lower_corner = [0] * ndims
        upper_corner = copy(cube_shape)

        # Set the channel
        lower_corner[spec_axis] = chan
        upper_corner[spec_axis] = chan

        box = rg.box(lower_corner, upper_corner)

        # Now make sliced image
        im_slice = ia.subimage(outputimage+"_channel_"+str(chan)+".image",
                               box)
        im_slice.done()


def get_slice_obj(slicearg):
    '''
    Given a string of the form "N:I:J", where N, I, J are integers,
    return a slice object: slice(N, I, J).

    http://stackoverflow.com/questions/680826/python-create-slice-object-from-string
    '''
    slice_ints = tuple([int(n) for n in slicearg.split(':')])
    return apply(slice, slice_ints)
