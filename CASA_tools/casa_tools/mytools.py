
'''
Define my own version of CASA tools that return python errors if they fail.
'''

from tasks import split, uvsub, concat, clean, rmtables, mstransform

from .graceful_error_catch import catch_fail


def mysplit(**kwargs):
    return catch_fail(split, **kwargs)


def myuvsub(**kwargs):
    return catch_fail(uvsub, **kwargs)


def myconcat(**kwargs):
    return catch_fail(concat, **kwargs)


def myclean(**kwargs):
    return catch_fail(clean, **kwargs)


def myrmtables(**kwargs):
    return catch_fail(rmtables, **kwargs)

def mymstransform(**kwargs):
    return catch_fail(mstransform, **kwargs)

