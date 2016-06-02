
from aplpy import FITSFigure
import os
import matplotlib.pyplot as mpl

'''
Show the moment 0 of the archival against the moment 0 from 14B-088
'''

fig = mpl.figure(figsize=(15, 7))

mom0_file = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits"

f1 = FITSFigure(mom0_file, figure=fig, subplot=[0.1,0.1,0.35,0.8])
# f1.set_tick_labels_font(size='x-small')
# f1.set_axis_labels_font(size='small')
f1.show_grayscale()

f1.hide_xaxis_label()
f1.hide_xtick_labels()

f1.hide_yaxis_label()
f1.hide_ytick_labels()

mom0_archival_file = "/media/eric/Data_3/M33/VLA_Data/AT0206/imaging/M33_206_b_c_HI.mom0.fits"

f2 = FITSFigure(mom0_archival_file, figure=fig, subplot=[0.5,0.1,0.35,0.8])
# f2.set_tick_labels_font(size='x-small')
# f2.set_axis_labels_font(size='small')
f2.show_grayscale()

f2.hide_xaxis_label()
f2.hide_xtick_labels()

f2.hide_yaxis_label()
f2.hide_ytick_labels()

fig.canvas.draw()
