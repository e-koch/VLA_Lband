
'''
Remove unneeded imaging products when only creating the dirty images.
'''

import os

from tasks import rmtables


def remove_products(name):
    '''
    For the dirty image, we don't need to keep the:

    * residual
    * flux
    * model

    If the image does not exist (i.e., no PSF due to flagging), remove all.
    '''

    removals = [".flux", ".residual", ".model"]

    # Check if the image exists
    if not os.path.exists(name + ".image"):
        rmtables(name + ".*")
    else:
        for remove in removals:
            rmtables(name + remove)
