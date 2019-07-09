
'''
Scripts for automatically setting clean parameters
'''

import numpy as np
from warnings import warn
import os

from taskinit import tb
from cleanhelper import cleanhelper


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

    def set_cellsize(vis, spw, sample_factor=6., baseline_percentile=95,
                     return_type="str"):

        syn_beam, prim_beam = \
            find_expected_beams(vis, spw,
                                baseline_percentile=baseline_percentile)

        # Round the cell size to some fraction, which becomes finer if it was
        # previously rounded to 0
        round_factor = 10
        while True:
            sel_cell_value = \
                round((syn_beam / sample_factor) * round_factor) / round_factor
            if sel_cell_value == 0:
                round_factor += 5
            else:
                break

        if return_type == "str":
            return str(sel_cell_value) + 'arcsec'
        else:
            return sel_cell_value

    def set_imagesize(vis, spw, source, sample_factor=6., pblevel=0.1,
                      max_size=15000, **kwargs):
        '''
        Set the image size for CLEAN to be a multiple of 2, 3, 5
        based on the maximum baseline in the MS.

        Parameters
        ----------
        '''

        if isinstance(max_size, (int, np.integer)):
            max_size = [max_size] * 2

        syn_beam, prim_beam = find_expected_beams(vis, spw)

        cellsize = set_cellsize(vis, spw, sample_factor=sample_factor,
                                return_type='value', **kwargs)

        if set_imagermode(vis, source) == "mosaic":
            mosaic_props = get_mosaic_info(vis, spw, sourceid=source,
                                           pblevel=pblevel)

            sel_imsize = [int(np.ceil(mosaic_props["Size_RA"] / cellsize)),
                          int(np.ceil(mosaic_props["Size_Dec"] / cellsize))]
        else:
            sel_imsize = [int(round(prim_beam / cellsize))] * 2

        # Check if this falls into the maximum allowed size. Otherwise, just
        # use the max.
        if sel_imsize[0] > max_size[0]:
            warn("Shape in first dimension exceeds maximum. Using maximum"
                 " given.")
            sel_imsize[0] = max_size[0]
        if sel_imsize[1] > max_size[1]:
            warn("Shape in second dimension exceeds maximum. Using maximum"
                 " given.")
            sel_imsize[1] = max_size[1]

        # The image size should be factorizable into some combo of
        # 2, 3, 5 and 7 to work with clean so:
        sel_imsize = [cleanhelper.getOptimumSize(size) for size in sel_imsize]

        # Return the rounded value nearest to the original image size chosen.
        return sel_imsize

    def find_expected_beams(vis, spw, baseline_percentile=95):
        '''
        Return the expected synthesized beam (approximately) and the primary
        beam size based on the baselines.

        Parameters
        ----------
        vis : str
            Name of MS.
        spw : int
            Which SPW in the MS to consider.
        baseline_percentile : int or float between 0 and 100
            The percentile of the longest baseline to estimate the synthesized
            beam with.

        Returns
        -------
        syn_beam : float
            Approximate Synthesized beam in arcseconds
        prim_beam : float
            Primary beam size in arcseconds.

        '''

        # Get percentile of max baseline and dish size
        bline_max = getBaselinePercentile(vis, baseline_percentile)

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
        centre_lambda = 299792458.0 / freq

        syn_beam = (centre_lambda / bline_max) * 180 / np.pi * 3600
        prim_beam = (centre_lambda / dish_min) * 180 / np.pi * 3600

        return syn_beam, prim_beam

    def getBaselinePercentile(msFile, percentile):
        """
        Based on getBaselineExtrema from analysisUtils
        """

        return np.percentile(getBaselines(msFile), percentile)

    def getBaselines(msFile):
        '''
        Return all baselines
        '''

        tb.open(msFile + '/ANTENNA')
        positions = np.transpose(tb.getcol('POSITION'))
        tb.close()

        all_lengths = []
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                length = au.computeBaselineLength(positions[i],
                                                  positions[j])
                if length != 0.0:
                    all_lengths.append(length)

        all_lengths = np.array(all_lengths)

        return all_lengths

    def get_mosaic_info(vis, spw, sourceid=None, intent='TARGET', pblevel=0.1):
        '''
        Return image size based on mosaic fields

        Parameters
        ----------
        vis : str
            MS Name
        spw : str or int
            If str, searches for an exact match with the names in the MS. If
            int, is the index of SPW in the MS.
        sourceid : str, optional
            The field names used will contain sourceid.
        intent : str, optional
            Use every field with the given intent (wildcards used by default).
        pblevel : float between 0 and 1
            PB level that defines the edges of the mosaic.
        '''

        mytb = au.createCasaTool(au.tbtool)

        # Check SPWs to make sure given choice is valid
        mytb.open(vis + '/SPECTRAL_WINDOW')

        spwNames = mytb.getcol('NAME')
        if isinstance(spw, str):
            match = False
            for spw_name in spwNames:
                if spw == spw_name:
                    match = True

            if not match:
                raise ValueError("The given SPW ({0}) is not in the MS SPW"
                                 " names ({1})".format(spw, spwNames))
        elif isinstance(spw, (int, np.integer)):
            try:
                spwNames[spw]
            except IndexError:
                raise IndexError("The given SPW index {0} is not in the range"
                                 " of SPWs in the MS ({1})."
                                 .format(spw, len(spwNames)))
        else:
            raise TypeError("spw must be a str or int.")

        refFreq = mytb.getcol("REF_FREQUENCY")[spw]
        lambdaMeters = au.c_mks / refFreq

        mytb.close()

        # Get field info
        mytb.open(vis + '/FIELD')

        delayDir = mytb.getcol('DELAY_DIR')
        ra = delayDir[0, :][0] * 12 / np.pi
        for i in range(len(ra)):
            if ra[i] < 0:
                ra[i] += 24
        ra *= 15
        dec = np.degrees(delayDir[1, :][0])

        # First choose fields by given sourceid
        if sourceid is not None:

            names = mytb.getcol('NAME')
            fields = mytb.getcol("SOURCE_ID")

            good_names = []
            good_fields = []
            for name, field in zip(names, fields):
                # Check if it has sourceid
                if name.find(sourceid) != -1:
                    good_names.append(name)
                    good_fields.append(field)
            names = good_names
            fields = good_fields

        # Then try choosing all fields based on the given intent.
        elif intent is not None:
            # Ensure a string
            intent = str(intent)

            mymsmd = au.createCasaTool(au.msmdtool)
            mymsmd.open(vis)
            intentsToSearch = '*' + intent + '*'
            fields = mymsmd.fieldsforintent(intentsToSearch)
            names = mymsmd.namesforfields(fields)
            mymsmd.close()
        # On or the other must be given
        else:
            raise ValueError("Either sourceid or intent must be given.")

        mytb.close()

        ra = ra[fields]
        dec = dec[fields]

        raAverageDegrees = np.mean(ra)
        decAverageDegrees = np.mean(dec)

        raRelativeArcsec = 3600 * (ra - raAverageDegrees) * \
            np.cos(np.deg2rad(decAverageDegrees))
        decRelativeArcsec = 3600 * (dec - decAverageDegrees)

        centralField = au.findNearestField(ra, dec,
                                           raAverageDegrees,
                                           decAverageDegrees)[0]

        # This next step is crucial, as it converts from the field number
        # determined from a subset list back to the full list.
        centralFieldName = names[centralField]
        centralField = fields[centralField]

        # Find which antenna have data
        mytb.open(vis)
        antennasWithData = np.sort(np.unique(mytb.getcol('ANTENNA1')))
        mytb.close()

        if antennasWithData.size == 0:
            raise Warning("No antennas with data found.")

        # Now we need the dish diameters
        mytb.open(vis + "/ANTENNA")
        # These are in m
        dish_diameters = \
            np.unique(mytb.getcol("DISH_DIAMETER")[antennasWithData])
        mytb.close()

        # Find maxradius
        maxradius = 0
        for diam in dish_diameters:
            arcsec = 0.5 * \
                au.primaryBeamArcsec(wavelength=lambdaMeters * 1000,
                                     diameter=diam, showEquation=False)
            radius = arcsec / 3600.0
            if radius > maxradius:
                maxradius = radius

        # Border about each point, down to the given pblevel
        border = 2 * maxradius * au.gaussianBeamOffset(pblevel) * 3600.

        size_ra = np.ptp(raRelativeArcsec) + 2 * border
        size_dec = np.ptp(decRelativeArcsec) + 2 * border

        mosaicInfo = {"Central_Field_ID": centralField,
                      "Central_Field_Name": centralFieldName,
                      "Center_RA": ra,
                      "Center_Dec": dec,
                      "Size_RA": size_ra,
                      "Size_Dec": size_dec}

        return mosaicInfo

    def append_to_cube(folder, prefix, suffix, num_imgs,
                       cube_name, chunk_size=250,
                       delete_chunk_cubes=True,
                       concat_kwargs={'relax': True, 'reorder': False,
                                      'overwrite': True}):
        '''
        Append single images into a cube. Must be continuous along
        spectral dimension.

        Individual images in the folder must be sequentially numbered.
        '''

        from os.path import join as osjoin
        import os

        try:
            import casatools
            ia = casatools.image()
        except ImportError:
            try:
                from taskinit import iatool
                ia = iatool()
            except ImportError:
                raise ImportError("Cannot import iatool.")

        imgs = [osjoin(folder, "{0}_{1}.{2}".format(prefix, chan, suffix))
                for chan in range(num_imgs)]

        # Make sure these all exist
        for img in imgs:
            if not os.path.exists(img):
                raise OSError("{} does not exist.".format(img))

        # Concatenate in chunks (Default of 250) b/c CASA doesn't like
        # having too many images open at once.

        num_chunks = (num_imgs // chunk_size) + 1

        chunk_cubes = []

        for i in range(num_chunks):

            start = chunk_size * i
            stop = min(chunk_size * (i + 1), num_imgs)

            imgs_chunk = imgs[start:stop]

            chunk_cube_name = "{0}_{1}".format(cube_name, i)
            chunk_cubes.append(chunk_cube_name)

            ia.imageconcat(outfile=chunk_cube_name,
                           infiles=imgs_chunk,
                           **concat_kwargs)
            ia.done()
            ia.close()

        # Concat the chunks together

        ia.imageconcat(outfile=cube_name,
                       infiles=chunk_cubes,
                       **concat_kwargs)
        ia.done()
        ia.close()

        if delete_chunk_cubes:
            for chunk_cube_name in chunk_cubes:
                os.system("rm -rf {}".format(chunk_cube_name))

except ImportError:
    warn("Could not import analysisUtils.")
