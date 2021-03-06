
from astropy.table import Table, Column
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as p
import numpy as np
from spectral_cube import SpectralCube

from analysis.paths import (iram_co21_data_path, fourteenB_HI_data_path,
                            paper1_figures_path)
from analysis.constants import moment0_name
from analysis.galaxy_params import gal

'''
Ignore "clouds" beyond 6.6 kpc to remove the edge effects
'''

tab = Table.read(iram_co21_data_path("m33.co21_new_props_clfind.fits"))

cloud_mask = \
    fits.open(iram_co21_data_path("m33.co21_new_assign_cprops.fits"))[0]

ico = fits.open(iram_co21_data_path("m33.ico.fits"))[0]

cube = SpectralCube.read(iram_co21_data_path("m33.co21_iram.fits"))
del cube._header[""]

radii = gal.radius(header=ico.header)

ax = p.subplot(111, projection=cube[0].wcs)

ax.imshow(ico.data.squeeze(), origin='lower', cmap='gray', vmin=0)
# p.colorbar()
# ax.contour((cloud_mask.data * (radii < 6500 * u.pc)).sum(0) > 0, colors='b')
ax.contour(radii < 6.5 * u.kpc, colors='b')
ax.scatter(tab["XPOS"], tab["YPOS"], transform=ax.get_transform('fk5'),
           edgecolor='r',
           facecolor='none')

p.draw()

p.savefig(paper1_figures_path("co21_moment0_cloud_centers_fullcat.pdf"))
p.savefig(paper1_figures_path("co21_moment0_cloud_centers_fullcat.png"))

# raw_input("Next plot?")
p.clf()

# Add a column for galactic radius
rgals = gal.radius(ra=(tab['XPOS']), dec=(tab['YPOS'])).to(u.kpc).value
colrgal = Column(name='RGAL_PC', data=(rgals))
tab.add_column(colrgal)

# Impose a 5 sigma cut
conds = np.logical_and(tab['MAXVAL'] / tab['NOISE'] > 5, tab["RGAL_PC"] < 6.5)
cleantab = tab[conds]

ax = p.subplot(111, projection=cube[0].wcs)

ax.imshow(ico.data.squeeze(), origin='lower', cmap='gray', vmin=0,
          vmax=np.nanpercentile(ico.data, 99))
# p.colorbar()
# ax.contour((cloud_mask.data * (radii < 6500 * u.pc)).sum(0) > 0, colors='b')
ax.contour(radii < 6.5 * u.kpc, colors='b')
ax.scatter(cleantab["XPOS"], cleantab["YPOS"],
           transform=ax.get_transform('fk5'), edgecolor='r',
           facecolor='none')

p.draw()

p.savefig(paper1_figures_path("co21_moment0_cloud_centers.pdf"))
p.savefig(paper1_figures_path("co21_moment0_cloud_centers.png"))

# raw_input("Next plot?")
p.clf()

# Save the clean sample
save = True
if save:
    cleantab.write(iram_co21_data_path("m33.co21_new_props_clfind_cleansample.fits"),
                   overwrite=True)

    # Remove and relabel the cloud mask.
    ctr = 0
    cloud_mask_data = cloud_mask.data.copy()
    for num, keep in enumerate(conds):
        # If considered a valid cloud, relabel
        if keep:
            cloud_mask_data[np.where(cloud_mask_data == num + 1)] = ctr
            ctr += 1
        # Else remove it.
        else:
            cloud_mask_data[np.where(cloud_mask_data == num + 1)] = 0

    # Save the new cloud mask.
    # This is in all caps. WCS doesn't like that.
    cloud_mask.header["CUNIT3"] = "km/s"
    new_cloud_mask = \
        fits.PrimaryHDU(cloud_mask_data,
                        header=WCS(cloud_mask.header).dropaxis(-1).to_header())
    new_cloud_mask.writeto(iram_co21_data_path("m33.co21_new_assign_cprops_cleansample.fits",
                                               no_check=True))


# Overplot the clean sample on the HI moment 0

mom0 = fits.open(fourteenB_HI_data_path(moment0_name))[0]

ax = p.subplot(111, projection=WCS(mom0.header))

ax.imshow(mom0.data, origin='lower', cmap='gray', vmin=0)
# p.colorbar()
# ax.contour((cloud_mask.data * (radii < 6500 * u.pc)).sum(0) > 0, colors='b')
hi_radii = gal.radius(header=mom0.header)
ax.contour(hi_radii < 6.5 * u.kpc, colors='b')
ax.scatter(cleantab["XPOS"], cleantab["YPOS"],
           transform=ax.get_transform('fk5'),
           edgecolor='r',
           facecolor='none')

p.draw()

p.savefig(paper1_figures_path("hi_moment0_cloud_centers.pdf"))
p.savefig(paper1_figures_path("hi_moment0_cloud_centers.png"))

# raw_input("Next plot?")
p.clf()
