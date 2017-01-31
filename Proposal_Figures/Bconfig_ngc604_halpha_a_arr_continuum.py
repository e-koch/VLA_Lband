
'''
Show preliminary continuum imaging from NGC 604 with Hubble images in
the Halpha.
'''

import os
import astropy.io.fits as fits
from aplpy import FITSFigure


halpha = fits.open("/mnt/MyRAID/M33/Hubble/ngc604/hst_05773_03_wfpc2_f656n_wf"
                   "/hst_05773_03_wfpc2_f656n_wf_sci.fits")[1]


lband_file = os.path.expanduser("~/bigdata/ekoch/VLA_tracks/reduction/"
                                "16B-242_10_16_16/16B-242_10_16_16_continuum/test_images/"
                                "NGC604_all_continuum.image.pbcor.fits")
lband_continuum = fits.open(lband_file)[0]

fig = FITSFigure(halpha)
fig.show_grayscale()

fig.show_contour(lband_continuum, levels=[3e-5, 7e-5, 1e-4, 5e-4, 7e-4, 1e-3,
                                          2e-3, 5e-3],
                 rasterize=True)

fig.hide_axis_labels()
fig.hide_tick_labels()

# Slice to immediate region about the cluster by-hand to avoid Hubble map
# edges showing.
raw_input("Adjust image extents.")

fig.save(os.path.expanduser("~/Dropbox/Various Plots/Proposals/NGC604_Halpha_A_config_continuum.pdf"))
fig.save(os.path.expanduser("~/Dropbox/Various Plots/Proposals/NGC604_Halpha_A_config_continuum.png"))

fig.close()
