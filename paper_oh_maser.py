
'''
Analysis and figures for research notes

Requires running OH_1665_narrowchannel_imaging
'''

from os.path import join as osjoin

from paths import c_path

# Masking and moment-making scripts

# Make velocity corrected cubes
execfile(osjoin(c_path, "Lines/OH_maser_luminosity.py"))
execfile(osjoin(c_path, "Lines/OH_maser_figure.py"))
