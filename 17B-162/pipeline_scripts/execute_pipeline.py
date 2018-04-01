
'''
Splits the ASDM into continuum and line portions, then runs the VLA pipeline
on each MS separately.

Copes the cont.dat file from this folder into the newly created MS folder.
'''

import os

script_path = None

execfile(os.path.join(script_path, "ms_split.py"))

orig_dir = os.getcwd()

os.chdir()
