
'''
Create PV-slices across the Northern "plume" region.
'''

from spectral_cube import SpectralCube, Projection
from astropy.io import fits
from astropy import units as u
import numpy as np
import pvextractor as pv
from regions import RectangleSkyRegion, write_ds9
from astropy.coordinates import Angle, SkyCoord
from aplpy import FITSFigure
import matplotlib.pyplot as plt


from paths import (fourteenB_HI_data_wGBT_path, allfigs_path,
                   fourteenB_wGBT_HI_file_dict)

rotsub_cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['RotSub_Cube'])

# Define perpendicular paths across the region

# in pixels
path1 = pv.Path([(580, 1320), (695, 1461)], width=15)
path2 = pv.Path([(680, 1300), (580, 1460)], width=15)


pvslice1 = pv.extract_pv_slice(rotsub_cube, path1)
pvslice2 = pv.extract_pv_slice(rotsub_cube, path2)

# Write out the paths as ds9 regions
x_cent1 = int((695 + 580) / 2.)
y_cent1 = int((1461 + 1320) / 2.)
x_cent2 = int((680 + 580) / 2.)
y_cent2 = int((1460 + 1300) / 2.)

pix_size = rotsub_cube.header['CDELT2'] * u.deg
height1 = np.sqrt((580 - 695)**2 + (1320 - 1461)**2) * pix_size
height2 = np.sqrt((580 - 680)**2 + (1300 - 1460)**2) * pix_size
width1 = 15 * pix_size
width2 = 15 * pix_size

angle1 = (((1461 - 1320) / (695. - 580.)) * u.rad).to(u.deg)
angle2 = (((1460 - 1300) / (580. - 680.)) * u.rad).to(u.deg)

cent1 = SkyCoord(rotsub_cube.spatial_coordinate_map[1][y_cent1, x_cent1],
                 rotsub_cube.spatial_coordinate_map[0][y_cent1, x_cent1], frame='icrs')
cent2 = SkyCoord(rotsub_cube.spatial_coordinate_map[1][y_cent2, x_cent2],
                 rotsub_cube.spatial_coordinate_map[0][y_cent2, x_cent2], frame='icrs')

rect_region1 = RectangleSkyRegion(center=cent1,
                                  height=Angle(height1),
                                  width=Angle(width1),
                                  angle=Angle(angle1))
rect_region2 = RectangleSkyRegion(center=cent2,
                                  height=Angle(height2),
                                  width=Angle(width2),
                                  angle=Angle(angle2))

write_ds9([rect_region1, rect_region2], allfigs_path("pvslices/nplume_pvslices.reg"))


# Now plot the pv-slices

fig1 = FITSFigure(pvslice1)
fig1.show_grayscale()
# I manually sliced these in the velocity direction before saving
p.tight_layout()
fig1.savefig(allfigs_path("pvslices/nplume_rotsub_pvslice1.pdf"))
fig1.savefig(allfigs_path("pvslices/nplume_rotsub_pvslice1.png"))
fig1.close()

fig2 = FITSFigure(pvslice2)
fig2.show_grayscale()
# I manually sliced these in the velocity direction before saving
p.tight_layout()
fig2.savefig(allfigs_path("pvslices/nplume_rotsub_pvslice2.pdf"))
fig2.savefig(allfigs_path("pvslices/nplume_rotsub_pvslice2.png"))
fig2.close()


# Now make a figure on the zeroth moment map of the two regions.

hdu = fits.open(fourteenB_wGBT_HI_file_dict['Moment0'])[0]

fig_mom0 = FITSFigure(hdu)
fig_mom0.show_grayscale()
# Note that I had to add metadata for the color/linestyle into the reg file
fig_mom0.show_regions(allfigs_path("pvslices/nplume_pvslices.reg"))
fig_mom0.savefig(allfigs_path("pvslices/nplume_pvslice_regions.pdf"))
fig_mom0.savefig(allfigs_path("pvslices/nplume_pvslice_regions.png"))
fig_mom0.close()

# Now make a projection so we can easily slice the zeroth moment
mom0 = Projection.from_hdu(hdu)

reg_slice = (slice(1250, 1550), slice(450, 750))
mom0_sliced = mom0[reg_slice]
mom0_sliced.quicklook()
mom0_sliced.FITSFigure.show_regions(allfigs_path("pvslices/nplume_pvslices.reg"))
mom0_sliced.FITSFigure.savefig(allfigs_path("pvslices/nplume_pvslice_regions_zoomed.pdf"))
mom0_sliced.FITSFigure.savefig(allfigs_path("pvslices/nplume_pvslice_regions_zoomed.png"))
mom0_sliced.FITSFigure.close()

# Now show a few rotsub channels in this zoom

fig = plt.figure()

fig1 = FITSFigure(rotsub_cube[1200][reg_slice].hdu, figure=fig, subplot=(2, 3, 1))
fig1.show_grayscale()
fig1.add_label(0.81, 0.92, "{} km/s".format(int(rotsub_cube.spectral_axis[1200].value / 1000.)), relative=True)
fig1.hide_axis_labels()
fig1.hide_tick_labels()

fig2 = FITSFigure(rotsub_cube[1150][reg_slice].hdu, figure=fig, subplot=(2, 3, 2))
fig2.show_grayscale()
fig2.add_label(0.81, 0.92, "{} km/s".format(int(rotsub_cube.spectral_axis[1150].value / 1000.)), relative=True)
fig2.hide_axis_labels()
fig2.hide_tick_labels()

fig3 = FITSFigure(rotsub_cube[1100][reg_slice].hdu, figure=fig, subplot=(2, 3, 3))
fig3.show_grayscale()
fig3.add_label(0.81, 0.92, "{} km/s".format(int(rotsub_cube.spectral_axis[1100].value / 1000.)), relative=True)
fig3.hide_axis_labels()
fig3.hide_tick_labels()

fig4 = FITSFigure(rotsub_cube[1050][reg_slice].hdu, figure=fig, subplot=(2, 3, 4))
fig4.show_grayscale()
fig4.add_label(0.81, 0.92, "{} km/s".format(int(rotsub_cube.spectral_axis[1050].value / 1000.)), relative=True)
fig4.hide_axis_labels()
fig4.hide_tick_labels()

fig5 = FITSFigure(rotsub_cube[1000][reg_slice].hdu, figure=fig, subplot=(2, 3, 5))
fig5.show_grayscale()
fig5.add_label(0.81, 0.92, "{} km/s".format(int(rotsub_cube.spectral_axis[1000].value / 1000.)), relative=True)
fig5.hide_axis_labels()
fig5.hide_tick_labels()

fig6 = FITSFigure(rotsub_cube[950][reg_slice].hdu, figure=fig, subplot=(2, 3, 6))
fig6.show_grayscale()
fig6.add_label(0.81, 0.92, "{} km/s".format(int(rotsub_cube.spectral_axis[950].value / 1000.)), relative=True)
fig6.hide_axis_labels()
fig6.hide_tick_labels()

p.tight_layout()

fig.savefig(allfigs_path("pvslices/nplume_rotsub_channels.png"))
fig.savefig(allfigs_path("pvslices/nplume_rotsub_channels.pdf"))
plt.close()
