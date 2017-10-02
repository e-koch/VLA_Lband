import numpy as np
import matplotlib.pyplot as p
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import average_beams
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import warnings
import matplotlib.animation as anim

from paths import fourteenB_HI_data_wGBT_path, allfigs_path, fourteenB_wGBT_HI_file_dict
from constants import hi_freq

'''
Channel plots of the rotation subtracted HI cube combined into a movie!

Borrowing code from @keflavich:
https://github.com/keflavich/paper_w51_evla/blob/master/plot_codes/h77a_layers.py
'''

cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['RotSube_Cube'])

# Begin channel map code here

# integrate over velocities to make channel maps of a set width
vstart = 0  # channels
vend = cube.shape[0]
vstep = 10
all_slabs = np.arange(vstart, vend + vstep, vstep, dtype=int)

# Define the average beam
try:
    beam = cube.beam
except AttributeError:
    beam = average_beams(cube.beams)

layers = \
    [cube[start:end].moment0().value *
     beam.jtok(hi_freq) / 1000. * u.km / u.s
     for start, end in
     ProgressBar(zip(all_slabs[:-1], all_slabs[1:]))]

# Scale all to the maximum

mx = np.max([np.nanmax(x).value for x in layers])

spec_axis = cube.spectral_axis.to(u.km / u.s).value

center_vels = [(spec_axis[start] + spec_axis[min(end, cube.shape[0] - 1)]) / 2. for start, end in
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

writer = anim.writers['ffmpeg'](fps=4)
ani.save(allfigs_path("m33_rotsub_movie.mp4"), writer=writer, dpi=300)
