
'''
Create a spectrum of the noise levels estimated when cleaning.
'''

import sys
import os
import numpy as np
from glob import glob
from astropy import units as u
from spectral_cube import (SpectralCube, OneDSpectrum)

path = sys.argv[1]
cube_name = sys.argv[2]

cube = SpectralCube.read(cube_name)

max_chan = cube.shape[0]

sigmas = []

stage1_level = 5.
stage2_level = 2.

for chan in range(max_chan):

    files = glob(os.path.join(path, "*channel_{}.*.npy".format(chan)))

    if len(files) == 0:
        raise ValueError("Missing npy files for channel {}".format(chan))

    is_stage2 = np.array(["stage2" in f for f in files])

    if is_stage2.any():
        chan_file = np.load(files[np.where(is_stage2)[0]]).item()
        try:
            sigma = (chan_file['summaryminor'][1][-2] / stage2_level)
            sigmas.append(sigma)
            continue

        except IndexError:
            files.remove(files[np.where(is_stage2)[0]])

    if not is_stage2.any():
        chan_file = np.load(files[0]).item()
        sigma = (chan_file['summaryminor'][1][-2] / stage1_level)
        sigmas.append(sigma)

spec = cube[:, 0, 0]

noise_spec = OneDSpectrum(sigmas, unit=u.Jy / u.beam, wcs=spec.wcs,
                          meta=spec.meta,
                          beams=cube.beams if hasattr(cube, 'beams') else None)

cube_path = os.path.dirname(cube_name)
cube_filename = os.path.split(cube_name)[-1]

noisespec_name = cube_filename[:-4] + "noise_spectrum.fits"

noise_spec.write(os.path.join(cube_path, noisespec_name), overwrite=True)
