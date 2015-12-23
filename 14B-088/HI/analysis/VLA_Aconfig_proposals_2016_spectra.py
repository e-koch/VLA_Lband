
'''
Extract interesting spectra from the full HI cube.
'''

import os
from spectral_cube import SpectralCube
import astropy.units as u
import matplotlib.pyplot as p
import seaborn as sn
sn.set_context('poster')
sn.set_style("ticks")

# path = "/srv/astro/erickoch/M33/"
path = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/"

cube = SpectralCube.read(os.path.join(path, "M33_14B-088_HI.clean.image.fits"),
                         mode='denywrite')

arec_cube = "/media/eric/Data_3/M33/Arecibo/14B-088_items/M33_14B-088_HI_model.fits"
# arec_cube = "/srv/astro/surveys/m33/hi/M33only.fits"

arecibo = SpectralCube.read(arec_cube, mode='denywrite')

# arec_spec = arecibo[:, 170, 132]

# I've imaged a portion of the Southern Arm w/o Arecibo.
vla_path = "/media/eric/MyRAID/M33/14B-088/HI/imaging/south_arm_800_1200.image.fits"
sarm_vla = SpectralCube.read(vla_path, mode='denywrite')

# Interesting positions by-eye
posns = [["South Arm 1", 1383, 1129, -210000, -75000],
         ["South Arm 2", 1427, 1099, -210000, -75000],
         ["South Arm 3", 1429, 1109, -210000, -75000],
         ["South Arm 4", 1413, 1106, -210000, -75000],
         ["South Arm 5", 1370, 1186, -210000, -75000],
         ["NGC 604 1", 1106, 1419, -280000, -150000],
         ["NGC 604 2", 1104, 1432, -280000, -150000],
         ["NGC 604 3", 1080, 1439, -280000, -150000],
         ["NGC 604 4", 1078, 1444, -280000, -150000]]

p.ioff()

for posn in posns:
    p.title(posn[0]+" at "+str(posn[1])+","+str(posn[2]))
    filename = \
        os.path.join("/home/eric/Dropbox/M33/",
                     "M33_14B-088_HI.clean.image." +
                     posn[0].replace(" ", "_").lower()+"_" +
                     str(posn[1])+"_"+str(posn[2])+".pdf")
    spec = cube[:, posn[2], posn[1]]
    spec_conv = spec.to(u.K, equivalencies=cube.beam.jtok_equiv(1420.40575177*u.MHz))
    # Temp fix for WCS not being kept after unit conversion
    spec_conv._wcs = spec.wcs
    spec_conv.quicklook(label="VLA + Arecibo", c='b')
    arec_spec = arecibo[:, -posn[2], -posn[1]]
    arec_spec_conv = arec_spec.to(u.K, equivalencies=arecibo.beam.jtok_equiv(1420.40575177*u.MHz))
    arec_spec_conv._wcs = arec_spec.wcs
    arec_spec_conv.quicklook(label="Arecibo", c='c')
    # p.plot(arec_spec_conv.spectral_axis[::-1], arec_spec_conv.value,
    #        drawstyle='steps-mid')
    if "South" in posn[0]:
        vla_spec = sarm_vla[:, posn[2], posn[1]]
        vla_spec_conv = vla_spec.to(u.K, equivalencies=sarm_vla.beam.jtok_equiv(1420.40575177*u.MHz))
        vla_spec_conv._wcs = vla_spec.wcs
        vla_spec_conv.quicklook(label="VLA", c='r')
    p.plot([posn[3], posn[4]], [0.0]*2, 'k--')
    p.xlim([posn[3], posn[4]])
    p.legend(loc='best')
    p.savefig(filename)
    p.show()
    p.clf()
