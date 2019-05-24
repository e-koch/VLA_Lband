
'''
Restore calibrated data from the pipeline products w/o hanning smoothing!

Does NOT use the standard folder setup. Assumes all relevant products are
in the current directory.


1) Start pipeline
2) Unpack flagversion -- Need to split out line SPWs so shouldn't have a
                         flagversions folder
3) Unpack cal tables
5) Apply last flag version to data (not necessarily 'Pipeline_Final')
4) Change directories in applycal text file. Feed into calstate in context
    - This ensures the tables and MS exist.
6) Apply cal tables.

** Tested on VLA pipeline in CASA 5.4.1. **
** Appears to work with the VLA pipeline products from 5.1.2 and above. **

'''

from __future__ import absolute_import

import os
import sys
import tempfile
import re
import tarfile

import pipeline
from pipeline.h.tasks import applycal as applycal_mod
import pipeline.infrastructure.vdp as vdp
import pipeline.infrastructure.tablereader as tablereader

from tasks import flagmanager


vis = sys.argv[-1]

__rethrow_casa_exceptions = True
context = h_init()

# context = pipeline.Pipeline().context
# inputs = pipeline.tasks.RestoreDataInputs(context, vis=vis)
# task = pipeline.tasks.RestoreData(inputs)
# results = task.execute(dry_run=False)
# results.accept(context)

try:

    # Untar the flags and cal tables

    flagname = "{}.flagversions.tgz".format(vis)
    if not os.path.exists(flagname):
        raise OSError("Cannot find flagversions.tgz: {}".format(flagname))
    with tarfile.open(flagname, 'r:gz') as tar:
        tar.extractall(path="")

    # Assume this is the name for now. Should be fine for all single
    # track pipeline runs
    tablename = "unknown.session_1.caltables.tgz"
    if not os.path.exists(tablename):
        raise OSError("Cannot find caltables.tgz: {}".format(tablename))
    with tarfile.open(tablename, 'r:gz') as tar:
        tar.extractall(path="")

    # Want to check if any custom flagging was made after
    # the `Pipeline_Final` version. If so, restore that version
    # instead.

    # in casa >5.4 (and quite possible earlier; I haven't checked),
    # flamanager returns a list of the versions
    out = flagmanager(vis, mode='list')

    # Last item is 'MS'. Don't need that
    items = out.keys()
    items.remove('MS')

    # Check the name of the last version
    last_flag_vs = out[max(items)]

    casalog.post("Restoring last flag version: "
                 "{}.".format(out[max(items)]['name']))

    flagmanager(vis, mode='restore',
                versionname=out[max(items)]['name'])

    # Register the MS with the pipeline session
    ms_reader = tablereader.ObservingRunReader

    observing_run = ms_reader.get_observing_run(vis)
    context.observing_run.add_measurement_set(observing_run.measurement_sets[0])

    # Update directories in applycal.txt
    applyfile = '{}.calapply.txt'.format(vis)

    if not os.path.exists(applyfile):
        raise OSError("Cannot find calapply file: {}".format(applyfile))

    unix_path = re.compile('((?:\\/[\\w\\.\\-]+)+)',
                           re.IGNORECASE | re.DOTALL)

    # Assume everything is in the same directory, so no osjoin in the return
    def repfn(matchobj):
        basename = os.path.basename(matchobj.group(0))
        return basename

    with open(applyfile, 'r') as f:
        out = unix_path.sub(repfn, f.read())

    with tempfile.NamedTemporaryFile() as tmpfile:
        tmpfile.write(out)
        tmpfile.flush()
        context.callibrary.import_state(tmpfile.name)

    # Run applycal
    container = vdp.InputsContainer(applycal_mod.Applycal, context)

    applycal_task = applycal_mod.Applycal(container)

    applycal_task.execute(dry_run=False)

    # could then run statwt. But going to disable that for the line
    # for now.
    # hifv_statwt()

finally:
    h_save()
