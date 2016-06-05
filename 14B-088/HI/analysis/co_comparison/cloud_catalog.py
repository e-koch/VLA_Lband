
from astropy.table import Table, Column
from astropy.io import fits
import astropy.units as u
import os
from galaxies import Galaxy
from astropy.wcs import WCS
import matplotlib.pyplot as p
import numpy as np
import aplpy
from spectral_cube import SpectralCube

'''
Ignore "clouds" beyond 6.6 kpc to remove the edge effects
'''

data_path = "/media/eric/Data_3/M33/co21/"

tab = Table.read(os.path.join(data_path, "m33.co21_new_props_clfind.fits"))

cloud_mask = fits.open(os.path.join(data_path, "m33.co21_new_assign_cprops.fits"))[0]

ico = fits.open(os.path.join(data_path, "m33.ico.fits"))[0]

cube = SpectralCube.read(os.path.join(data_path, "m33.co21_iram.fits"))

gal = Galaxy("M33")
radii = gal.radius(header=ico.header)

ax = p.subplot(111, projection=cube[0].wcs)

ax.imshow(ico.data.squeeze(), origin='lower', cmap='gray', vmin=0)
# p.colorbar()
# ax.contour((cloud_mask.data * (radii < 6500 * u.pc)).sum(0) > 0, colors='b')
ax.contour(radii < 6.5 * u.kpc, colors='b')
ax.scatter(tab["XPOS"], tab["YPOS"], transform=ax.get_transform('fk5'), edgecolor='r',
           facecolor='none')

p.draw()
raw_input("Next plot?")
p.clf()

# Add a column for galactic radius
rgals = gal.radius(ra=(tab['XPOS']), dec=(tab['YPOS'])).to(u.kpc).value
colrgal=Column(name='RGAL_PC', data=(rgals))
tab.add_column(colrgal)

# Impose a 5 sigma cut
conds = np.logical_and(tab['MAXVAL'] / tab['NOISE'] > 5, tab["RGAL_PC"] < 6.5)
cleantab = tab[conds]

ax = p.subplot(111, projection=cube[0].wcs)

ax.imshow(ico.data.squeeze(), origin='lower', cmap='gray', vmin=0)
# p.colorbar()
# ax.contour((cloud_mask.data * (radii < 6500 * u.pc)).sum(0) > 0, colors='b')
ax.contour(radii < 6.5 * u.kpc, colors='b')
ax.scatter(cleantab["XPOS"], cleantab["YPOS"], transform=ax.get_transform('fk5'), edgecolor='r',
           facecolor='none')

# Save the clean sample
save = True
if save:
    cleantab.write(os.path.join(data_path, "m33.co21_new_props_clfind_cleansample.fits"))

p.draw()
raw_input("Next plot?")
p.clf()

# Overplot the clean sample on the HI moment 0

mom0_file = "/media/eric/MyRAID/M33/14B-088/HI/full_imaging/M33_14B-088_HI.clean.image.pbcov_gt_0.3_masked.mom0.fits"
mom0 = fits.open(mom0_file)[0]

ax = p.subplot(111, projection=WCS(mom0.header))

ax.imshow(mom0.data, origin='lower', cmap='gray', vmin=0)
# p.colorbar()
# ax.contour((cloud_mask.data * (radii < 6500 * u.pc)).sum(0) > 0, colors='b')
hi_radii = gal.radius(header=mom0.header)
ax.contour(hi_radii < 6.5 * u.kpc, colors='b')
ax.scatter(cleantab["XPOS"], cleantab["YPOS"], transform=ax.get_transform('fk5'),
           edgecolor='r',
           facecolor='none')
