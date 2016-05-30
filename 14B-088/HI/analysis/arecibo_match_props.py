
from spectral_cube import SpectralCube
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.convolution import convolve_fft
from astropy.utils.console import ProgressBar
import astropy.units as u
import numpy as np
import matplotlib.pyplot as p

'''
Compare properties of the VLA + Arecibo and VLA + Archival VLA + Arecibo to
the Arecibo only

I'm using the Arecibo cube upsampled onto the VLA grid, which was used as the
model during cleaning.
'''

verbose = False
which_cube = "new"

if which_cube not in ["new", "combined", "archival"]:
    raise ValueError("Bad which_cube: {}".format(which_cube))

if which_cube == "new":
    cube = SpectralCube.read("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits")
    arecibo = SpectralCube.read("/media/eric/Data_3/M33/Arecibo/14B-088_items/M33_14B-088_HI_model.fits")
elif which_cube == "combined":
    cube = SpectralCube.read("/media/eric/MyRAID/M33/14B-088/HI/combined_HI/M33_14B-088_AT0206_HI.clean.image.pbcov_gt_0.3_masked.fits")
    arecibo = SpectralCube.read("/media/eric/Data_3/M33/Arecibo/AT0206_items/M33_model.fits")
else:
    cube = SpectralCube.read("/media/eric/Data_3/M33/VLA_Data/AT0206/imaging/M33_206_b_c_HI.fits")
    arecibo = SpectralCube.read("/media/eric/Data_3/M33/Arecibo/AT0206_items/M33_model.fits")

# Now cut Arecibo to match the shapes.

arecibo = arecibo.spectral_slab(*cube.spectral_extrema)
arecibo = arecibo.subcube(ylo=cube.latitude_extrema[1],
                          yhi=cube.latitude_extrema[0],
                          xlo=cube.longitude_extrema[0],
                          xhi=cube.longitude_extrema[1])

arec_sums = np.empty((cube.shape[0], ))
vla_sums = np.empty((cube.shape[0], ))
vla_conv_sums = np.empty((cube.shape[0], ))

for i in ProgressBar(xrange(cube.shape[0])):
    # Sum over a channel
    ar = arecibo[i].sum() * \
        (1 / arecibo.beam.sr.to(u.deg ** 2)) * \
        (arecibo.header["CDELT2"] * u.deg)**2

    vla = np.nansum(cube[i]) * \
        (1 / cube.beam.sr.to(u.deg ** 2)) * \
        (cube.header["CDELT2"] * u.deg)**2

    if verbose:
        print("Total sum at native resolution.")
        print("Arecibo: {0}, VLA + Arecibo: {1}".format(ar.value, vla.value))
        # print("Combined VLA + Arecibo: {}".format(vla_combined.value))

    arec_sums[i] = ar.value
    vla_sums[i] = vla.value
    # vla_combined_sums[i] = vla_combined.value

    # Now convolve VLA to match Arecibo
    # This should match the previous sum

    pixscale = proj_plane_pixel_scales(arecibo.wcs)[0]
    beam_kern = arecibo.beam.as_kernel(pixscale)

    # conv_channel = convolve_fft(cube[i].value, beam_kern)

    # vla_conv = np.nansum(conv_channel) * \
    #     (1 / arecibo.beam.sr.to(u.deg ** 2)) * \
    #     (arecibo.header["CDELT2"] * u.deg)**2

    # if verbose:
    #     print("Total sum at common resolution.")
    #     print("Arecibo: {0}, VLA + Arecibo: {1}".format(ar.value, vla_conv.value))

    #     print("VLA + Arecibo should be scaled by:"
    #           " {}".format(vla_conv.value / ar.value))

    # vla_conv_sums[i] = vla_conv

# Now let's plot the channel sums
chans = np.arange(cube.shape[0])
p.plot(chans, arec_sums, 'b-', label='Arecibo', drawstyle='steps-mid')
p.plot(chans, vla_sums, 'g--', label='Arecibo + VLA', drawstyle='steps-mid')
# p.plot(chans, vla_conv_sums, 'k-.', label='Convolved VLA')
p.legend()
p.xlabel("Channel")
p.ylabel("Total Intensity (Jy)")
p.show()

# raw_input("Next plot?")

# p.plot(chans, arec_sums, 'b-', label='Arecibo')
# p.plot(chans, vla_sums / 1.451, 'g--', label='Arecibo + VLA', drawstyle='steps-mid')
# # p.plot(chans, vla_conv_sums, 'k-.', label='Convolved VLA', drawstyle='steps-mid')
# p.legend()
# p.xlabel("Channel")
# p.ylabel("Total Intensity (Jy)")
# p.show()

# Up to 973, where Arecibo seems to have more emission coming from the HI cloud
# that is not as well recovered in the VLA cube.
if which_cube == "new":
    slicer = slice(0, 973)
else:
    slicer = slice(None)
# Remove any zeros in the Arecibo
arec_zeros = np.nonzero(arec_sums)

avg_scale_factor = np.nanmean(vla_sums[arec_zeros][slicer] / arec_sums[arec_zeros][slicer])
std_scale_factor = np.nanstd(vla_sums[arec_zeros][slicer] / arec_sums[arec_zeros][slicer])
print("Scale Factor is: {0}+/-{1}".format(avg_scale_factor, std_scale_factor))
# New cube : 1.451 +/- 0.036
# Combined cube : 1.16695609064+/-0.239566349339
# Archival cube : 1.71537328973+/-1.38680301954  Eew. Largely driven by a
# couple of weird channels, which honestly should probably be reimaged.
