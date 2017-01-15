
import os
import matplotlib.pyplot as plt
import numpy as np

from plotting_styles import onecolumn_figure, default_figure
from paths import paper1_figures_path

'''
Make a UV plot of the 1000th HI channel.
'''

uvw = np.load("/mnt/MyRAID/M33/VLA/14B-088/HI/"
              "14B-088_HI_LSRK.ms.contsub_channel_1000.uvw.npy")

onecolumn_figure()
fig = plt.figure()
ax = fig.add_subplot(111)  # , rasterized=True)
# plt.hexbin(uvw[0], uvw[1], bins='log', cmap='afmhot_r')
ax.scatter(uvw[0], uvw[1], s=0.1, color='k', rasterized=True)
plt.xlabel("U (m)")
plt.ylabel("V (m)")
plt.xlim([-3200, 3500])
plt.ylim([-3200, 3200])
plt.grid()
plt.tight_layout()

plt.savefig(paper1_figures_path("m33_hi_uv_plane_chan1000.pdf"))
plt.savefig(paper1_figures_path("m33_hi_uv_plane_chan1000.png"))

plt.close()

default_figure()
