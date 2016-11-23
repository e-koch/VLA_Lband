import numpy as np
import os
import matplotlib.pyplot as p
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import average_beams
from astropy import coordinates
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy import log
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import warnings
import matplotlib.animation as anim

from analysis.paths import fourteenB_HI_data_path, paper1_figures_path

'''
Channel plots of the rotation subtracted HI cube combined into a movie!

Borrowing code from @keflavich: https://github.com/keflavich/paper_w51_evla/blob/master/plot_codes/h77a_layers.py
'''

cube_file = fourteenB_HI_data_path("M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.rotsub.fits")

cube = SpectralCube.read(cube_file)

# Begin channel map code here

# integrate over velocities to make channel maps of a set width
vstart = 308  # channels
vend = 788
vstep = 10
all_slabs = np.arange(vstart, vend + vstep, vstep, dtype=int)

# Define the average beam
beam = average_beams(cube.beams)

layers = \
    [cube[start:end].moment0().value *
     beam.jtok(1.414 * u.GHz) / 1000. * u.km / u.s
     for start, end in zip(all_slabs[:-1], all_slabs[1:])]

# Scale all to the maximum
mx = np.max([np.nanmax(x).value for x in layers])

spec_axis = cube.spectral_axis.to(u.km / u.s).value

center_vels = [(spec_axis[start] + spec_axis[end]) / 2. for start, end in
               zip(all_slabs[:-1], all_slabs[1:])]

pb = ProgressBar(len(center_vels))

fig = p.figure()
ax = fig.add_subplot(111)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

p.tight_layout()


def updater(i):

    pb.update()

    layer = layers[i]

    im = ax.imshow(layer.value, origin='lower',
                   norm=ImageNormalize(vmin=-0.001,
                                       vmax=mx,
                                       stretch=AsinhStretch()),
                   cmap=p.cm.gray_r)
    # ax.annotate("${0:.0f} km/s$".format(center_vels[i]),
    #             (0.53, 0.9),
    #             xycoords='axes fraction', color='k',
    #             fontsize=15.5)


ani = anim.FuncAnimation(fig, updater, range(len(center_vels)))
# p.show()

p.ioff()
writer = anim.writers['ffmpeg'](fps=3)
ani.save(paper1_figures_path("m33_rotsub_movie.mp4"), writer=writer, dpi=300)
p.ion()
