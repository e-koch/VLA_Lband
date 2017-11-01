
'''
Adventures in clustering the HI spectra where CO is detected.
'''

from spectral_cube import SpectralCube, Projection
import numpy as np
import pylab as pl
import time
from astropy.io import fits
from sklearn.cluster import AgglomerativeClustering

from paths import fourteenB_wGBT_HI_file_dict, iram_co21_14B088_data_path
from plotting_styles import default_figure, onecolumn_figure, twocolumn_figure

# cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['PeakSub_Cube'])
cube = SpectralCube.read(fourteenB_wGBT_HI_file_dict['Cube'])
co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI.fits"))
# co_cube = SpectralCube.read(iram_co21_14B088_data_path("m33.co21_iram.14B-088_HI_feather.peakvels_corrected.fits"))

co_mask = fits.open(iram_co21_14B088_data_path(
    "m33.co21_iram.14B-088_HI_source_mask.fits"))[0].data
co_mask = co_mask.sum(0) > 3
yposns, xposns = np.where(co_mask)

spectra = np.vstack([cube[:, y, x].value for y, x in zip(yposns, xposns)])

# Be lazy. Use default parameters.
model = AgglomerativeClustering(linkage='average',
                                connectivity=None,
                                n_clusters=2)
t0 = time.time()
model.fit(spectra)
elapsed_time = time.time() - t0
print(elapsed_time)

class0_spec = spectra[model.labels_ == 0]
class1_spec = spectra[model.labels_ == 1]

onecolumn_figure()
pl.plot(class0_spec.mean(0), label='Class 0')
pl.plot(class1_spec.mean(0), label='Class 1')

# Where are these spatially located?
peaktemps = Projection.from_hdu(fits.open(fourteenB_wGBT_HI_file_dict['PeakTemp'])[0])

twocolumn_figure()
peaktemps.quicklook(use_aplpy=False)

pl.plot(xposns[model.labels_ == 1],
        yposns[model.labels_ == 1], 'ro', alpha=0.5)
pl.plot(xposns[model.labels_ == 1],
        yposns[model.labels_ == 1], 'bD', alpha=0.7)

default_figure()
