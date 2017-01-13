
import numpy as np
from taskinit import mstool


def extract_ms_data(vis, columns, field="", **kwargs):
    '''
    Extract data points from an MS.
    '''

    myms = mstool()

    myms.open(vis)
    myms.selectinit(0)
    selection_dict = dict(field=field)
    selection_dict.update(kwargs)
    assert myms.msselect(selection_dict), "Data selection has failed"

    if not isinstance(columns, list):
        columns = list(columns)

    datadict = myms.getdata(columns)
    myms.close()

    for column in columns:
        np.save("{0}_{1}.npy".format(vis.rstrip(".ms"), column.lower()),
                datadict[column])
