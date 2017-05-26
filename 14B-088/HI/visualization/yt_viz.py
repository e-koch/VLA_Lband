
'''
Image rendering with yt's Scene functionality.
'''

from spectral_cube import SpectralCube
from astropy.io import fits
import yt
import numpy as np

cube = SpectralCube.read("M33_14B-088_HI.clean.image.pbcov_gt_0.5_masked.fits")
mask = fits.open("M33_14B-088_HI.clean.image.pbcov_gt_0.5_masked_source_mask.fits")[0]

cube = cube.with_mask(mask.data > 0)

cube_max = cube.max()

ytcube = cube[100:-100].to_yt()
ytcube.dataset.periodicity = (True, True, True)

sc = yt.create_scene(ytcube.dataset, field='flux')
cam = sc.camera
source = sc[0]
source.log_field = True
source.set_use_ghost_zones(True)

cam.resolution = [2048, 2048]

# Start with the disk face on.
cam.rotate(np.pi / 2.)

bounds = (np.log10(0.004), np.log10(cube_max.value))
# bounds = (0.001, cube.max().value)
tf = yt.ColorTransferFunction(bounds)

# def linramp(vals, minval, maxval):
#     return (vals - vals.min()) / (vals.max() - vals.min())


# tf.map_to_colormap(bounds[0], bounds[1],
#                    colormap='arbre',
#                    scale_func=linramp)
nLayer = 10
# tf.add_layers(nLayer, colormap='gist_rainbow',
#               alpha=np.linspace(0, 1, nLayer))
tf.add_layers(nLayer, w=0.001, colormap='gist_rainbow',
              alpha=np.logspace(-1.0, -0.2, nLayer))

source.tfh.tf = tf
source.tfh.set_log(False)
source.tfh.bounds = bounds


# save an image at the starting position
frame = 0
sc.save('yt_outputs/camera_movement_%04i.png' % frame,
        sigma_clip=3)
frame += 1

# for _ in cam.iter_zoom(1.5, 2):
#     sc.render()
#     sc.save('yt_outputs/camera_movement_%04i.png' % frame,
#             sigma_clip=1)
#     frame += 1

for _ in cam.iter_rotate(1.9 * np.pi, 20):
    sc.render()
    sc.save('yt_outputs/camera_movement_%04i.png' % frame,
            sigma_clip=1)
    frame += 1
