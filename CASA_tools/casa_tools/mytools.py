
'''
Define my own version of CASA tools that return python errors if they fail.
'''

from tasks import split, uvsub, concat, clean, rmtables, mstransform

from taskinit import ia

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
