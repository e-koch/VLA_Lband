
'''
Scripts for automatically setting clean parameters
'''

import numpy as np
from warnings import warn
import os

from taskinit import tb


def set_imagermode(vis, source):

    tb.open(os.path.join(vis, 'FIELD'))
    names = tb.getcol('NAME')
    tb.close()

    moscount = 0

    for name in names:
        chsrc = name.find(source)

        if chsrc != -1:
            moscount = moscount + 1

    if moscount > 1:
        imagermode = "mosaic"
    else:
        imagermode = "csclean"

    return imagermode


def has_field(vis, source):
    '''
    Check if source is contained in at least one of the field names.
    '''

    tb.open(os.path.join(vis, 'FIELD'))
    names = tb.getcol('NAME')
    tb.close()

    moscount = 0

    for name in names:
        chsrc = name.find(source)

        if chsrc != -1:
            moscount = moscount + 1

    if moscount == 0:
        return False

    return True

try:
    import analysisUtils as au

    def set_cellsize(vis, spw, sample_factor=6.):

        syn_beam, prim_beam = find_expected_beams(vis, spw)

        # Setting CLEANing parameters
        sel_cell = str(round((syn_beam / sample_factor) * 2) / 2) + \
            'arcsec'

        return sel_cell

    def set_imagesize(vis, spw, sample_factor=6.):
        '''
        Set the image size for CLEAN to be a multiple of 2, 3, 5
        based on the maximum baseline in the MS.
        '''

        syn_beam, prim_beam = find_expected_beams(vis, spw)

        sel_imsize = int(round(prim_beam / (syn_beam / sample_factor)))

        # Increase the sel_imsize by a couple of beams
        # to be sure
        dx = int(round(syn_beam / prim_beam * sel_imsize))
        sel_imsize = sel_imsize + 1 * dx

        # The image size should be a multiplier of
        # 2, 3 and 5 to work well with clean so:

        rounded_sizes = []
        for base in [2, 5, 10]:
            rounded_sizes.append(round_to_base(sel_imsize, base=base))
        rounded_sizes = np.array(rounded_sizes)
        nearest_idx = np.abs(rounded_sizes - sel_imsize).argmin()

        # Return the rounded value nearest to the original image size chosen.
        return rounded_sizes[nearest_idx]

    def round_to_base(x, base=5):
        return int(base * round(float(x) / base))

    def find_expected_beams(vis, spw):
        '''
        Return the expected synthesized beam (approximately) and the primary
        beam size based on the baselines.

        Parameters
        ----------
        vis : str
            Name of MS.
        spw : int
            Which SPW in the MS to consider.

        Returns
        -------
        syn_beam : float
            Approximate Synthesized beam in arcseconds
        prim_beam : float
            Primary beam size in arcseconds.

        '''

        # Get max baseline and dish size
        bline_max = au.getBaselineExtrema(vis)[0]

        tb.open(os.path.join(vis, 'ANTENNA'))
        dishs = tb.getcol('DISH_DIAMETER')
        dish_min = min(dishs)
        tb.close()

        tb.open(os.path.join(vis, 'SPECTRAL_WINDOW'))
        ref_freqs = tb.getcol('REF_FREQUENCY')

        try:
            freq = ref_freqs[spw]
        except IndexError:
            raise IndexError("Given SPW ({0}) is not within the range of SPWs"
                             "found ({1})".format(spw, len(ref_freqs)))

        # Find the beam
        # XXX
        # When astropy is easier to install (CASA 4.7??), switch to using the
        # defined constants.
        centre_lambda = 299792458.0 / (freq)
        # min_lambda = 299792458.0 / (min(ref_freqs))
        syn_beam = (centre_lambda / bline_max) * 180 / np.pi * 3600
        prim_beam = (centre_lambda / dish_min) * 180 / np.pi * 3600

        return syn_beam, prim_beam

except ImportError:
    warn("Could not import analysisUtils.")
