
'''
Extract interesting spectra from the full HI cube.
'''

import os
from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as p
import seaborn as sn
sn.set_context('poster')
sn.set_style("dark")

path = "/srv/astro/erickoch/M33/"

cube = SpectralCube.read(os.path.join(path, "M33_14B-088_HI.clean.image.fits"))

# Interesting positions by-eye
posns = [["South Arm 1", 1383, 1129],
         ["South Arm 2", 1427, 1099],
         ["South Arm 3", 1429, 1109],
         ["South Arm 4", 1413, 1106],
         ["South Arm 5", 1370, 1186],
         ["NGC 604 1", 1106, 1419],
         ["NGC 604 2", 1104, 1432],
         ["NGC 604 3", 1080, 1439],
         ["NGC 604 4", 1078, 1444]]

for posn in posns:
    p.title(posn[0]+" at "+str(posn[1])+","+str(posn[2]))
    filename = \
        os.path.join(path,
                     "M33_14B-088_HI.clean.image." +
                     posn[0].replace(" ", "_").lower()+"_" +
                     str(posn[1])+"_"+str(posn[2]))
    cube[:, posn[1], posn[2]].quicklook(filename=filename)
