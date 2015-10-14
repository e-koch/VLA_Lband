
'''
Split MS into single spectral channels
'''

import os

from taskinit import tb

from .mytools import mymstransform


def split_by_channel(vis, nchan=-1, start=1, spw='0',
                     restfreq='1420.40575177MHz', **kwargs):
    '''
    Splits a MS by its spectral channels, according to the given
    '''

    # Get the total number of channels
    tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
    nchan_list = tb.getcol('NUM_CHAN')
    tb.close()

    # Pick out the spw specified
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

    for chan in channel_list:
        channel_vis = vis.rstrip(".ms")+"_channel_"+str(chan)+".ms"
        mymstransform(vis=vis, outputvis=channel_vis, spw=spw,
                      regridms=True, width=1, nchan=1, start=chan,
                      restfreq=restfreq, **kwargs)


def get_slice_obj(slicearg):
    '''
    Given a string of the form "N:I:J", where N, I, J are integers,
    return a slice object: slice(N, I, J).

    http://stackoverflow.com/questions/680826/python-create-slice-object-from-string
    '''
    slice_ints = tuple([int(n) for n in slicearg.split(':')])
    return apply(slice, slice_ints)
