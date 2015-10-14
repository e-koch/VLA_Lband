
'''
Split MS into single spectral channels
'''

import os

from taskinit import tb

from .mytools import mymstransform


def split_by_channel(vis, nchan=-1, start=1, spw=0,
                     restfreq='1420.40575177MHz'):
    '''
    Splits a MS by its spectral channels, according to the given
    '''

    # Get the total number of channels
    tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
    nchan_list = tb.getcol('NUM_CHAN')
    tb.close()

    # Pick out the spw specified
    total_nchan = nchan_list[spw]

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
                      restfreq=restfreq)
