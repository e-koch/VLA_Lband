
'''
Define my own version of CASA tools that return python errors if they fail.
'''

from tasks import split, uvsub, concat, clean, rmtables, mstransform, \
    exportfits

from taskinit import ia

from astropy.io import fits
import astropy.units as u

from .graceful_error_catch import catch_fail


def mysplit(**kwargs):
    return catch_fail(split, **kwargs)


def myuvsub(**kwargs):
    return catch_fail(uvsub, **kwargs)


def myconcat(**kwargs):
    return catch_fail(concat, **kwargs)


def myrmtables(**kwargs):
    return catch_fail(rmtables, **kwargs)


def mymstransform(**kwargs):
    return catch_fail(mstransform, **kwargs)


def myclean(**kwargs):

    if "mask" in kwargs.keys():
        # Since masks can be given in other forms, try to open it, and if it
        # fails, assume it just isn't an image.
        try:
            # Check if there's anything in the mask before cleaning
            ia.open(kwargs["mask"])
            stats_dict = ia.statistics()
            ia.close()
            # If there's nothing there, max == min
            max_val = stats_dict["max"]
            if max_val == 0:
                Warning("The mask image contains no regions to clean in."
                        " Exiting.")
                return None
        except RuntimeError:
            pass

    return catch_fail(clean, **kwargs)


def myexportfits(common_beam=True, **kwargs):
    '''
    Version of exportfits that returns a Python error when it fails.
    Also attachs a common beam to the fits header.
    '''

    catch_fail(exportfits, **kwargs)

    # The above throws an error if exportfits fails.
    # Now open the FITS file and attach a beam.

    ia.open(kwargs['imagename'])

    com_beam = ia.commonbeam()

    ia.close()

    bmaj = com_beam['major']['value'] * \
        u.Unit(com_beam['major']['unit']).to(u.deg)
    bmin = com_beam['minor']['value'] * \
        u.Unit(com_beam['minor']['unit']).to(u.deg)
    bpa = com_beam['pa']['value'] * u.Unit(com_beam['pa']['unit']).to(u.deg)

    filename = kwargs['fitsimage']

    output_fits = fits.open(filename, mode='update')

    output_fits[0].header.update({"BMAJ": bmaj,
                                  "BMIN": bmin,
                                  "BPA": bpa})

    output_fits.flush()
    output_fits.close()
