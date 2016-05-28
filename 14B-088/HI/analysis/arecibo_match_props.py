
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

cube = SpectralCube.read("/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.fits")
comb_cube = SpectralCube.read("/media/eric/MyRAID/M33/14B-088/HI/combined_HI/M33_14B-088_AT0206_HI.clean.image.pbcov_gt_0.3_masked.fits")
arecibo = SpectralCube.read("/media/eric/Data_3/M33/Arecibo/14B-088_items/M33_14B-088_HI_model.fits")

# Now cut Arecibo to match the shapes.

arecibo = arecibo.spectral_slab(*cube.spectral_extrema)
arecibo = arecibo.subcube(ylo=cube.latitude_extrema[1],
                          yhi=cube.latitude_extrema[0],
                          xlo=cube.longitude_extrema[0],
                          xhi=cube.longitude_extrema[1])
comb_cube = comb_cube.subcube(ylo=cube.latitude_extrema[0],
                              yhi=cube.latitude_extrema[1],
                              xlo=cube.longitude_extrema[1],
                              xhi=cube.longitude_extrema[0])

arec_sums = np.empty((cube.shape[0], ))
vla_sums = np.empty((cube.shape[0], ))
vla_combined_sums = np.empty((cube.shape[0], ))
vla_conv_sums = np.empty((cube.shape[0], ))

for i in ProgressBar(xrange(cube.shape[0])):
    # Sum over a channel
    ar = arecibo[i].sum() * \
        (1 / arecibo.beam.sr.to(u.deg ** 2)) * \
        (arecibo.header["CDELT2"] * u.deg)**2

    vla = np.nansum(cube[i]) * \
        (1 / cube.beam.sr.to(u.deg ** 2)) * \
        (cube.header["CDELT2"] * u.deg)**2

    vla_combined = np.nansum(comb_cube[i]) * \
        (1 / comb_cube.beam.sr.to(u.deg ** 2)) * \
        (comb_cube.header["CDELT2"] * u.deg)**2

    if verbose:
        print("Total sum at native resolution.")
        print("Arecibo: {0}, VLA + Arecibo: {1}".format(ar.value, vla.value))
        print("Combined VLA + Arecibo: {}".format(vla_combined.value))

    arec_sums[i] = ar.value
    vla_sums[i] = vla.value
    vla_combined_sums[i] = vla_combined.value

    # Now convolve VLA to match Arecibo
    # This should match the previous sum

    pixscale = proj_plane_pixel_scales(arecibo.wcs)[0]
    beam_kern = arecibo.beam.as_kernel(pixscale)

    conv_channel = convolve_fft(cube[i].value, beam_kern)

    vla_conv = np.nansum(conv_channel) * \
        (1 / arecibo.beam.sr.to(u.deg ** 2)) * \
        (arecibo.header["CDELT2"] * u.deg)**2

    if verbose:
        print("Total sum at common resolution.")
        print("Arecibo: {0}, VLA + Arecibo: {1}".format(ar.value, vla.value))

        print("VLA + Arecibo should be scaled by:"
              " {}".format(vla.value / ar.value))

    vla_conv_sums[i] = vla_conv

# Now let's plot the channel sums
chans = np.arange(cube.shape[0])
p.plot(chans, arec_sums, 'b-', label='Arecibo')
p.plot(chans, vla_sums, 'g--', label='Arecibo + VLA')
p.plot(chans, vla_sums, 'c.', label='Arecibo + VLA + Archival VLA')
p.plot(chans, vla_conv_sums, 'k-.', label='Convolved VLA')
p.legend()
p.xlabel("Channel")
p.ylabel("Total Intensity (Jy)")
p.show()
