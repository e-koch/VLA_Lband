
# Convert the CASA cube into FITS

import os
from spectral_cube import SpectralCube
import astropy.units as u

from casa_tools import myexportfits

path = "/srv/astro/erickoch/M33/"

myexportfits(imagename=os.path.join(path, "M33_14B-088_HI.clean.image"),
             fitsimage=os.path.join(path, "M33_14B-088_HI.clean.image.fits"),
             velocity=True, dropstokes=True, history=False)

cube = SpectralCube.read(os.path.join(path, "M33_14B-088_HI.clean.image.fits"),
                         mode='update')

converted_cube = cube.to(u.K, equivalencies=cube.beam.jtok_equiv(1.42040575177*u.GHz))

converted_cube.write(os.path.join(path, "M33_14B-088_HI.clean.image.fits"),
                     overwrite=True)
