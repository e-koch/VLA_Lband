
'''
Investigate simple spatial distributions of discrete structures to
determine if the qualitative properties of the CO(2-1) can be reproduced
'''


import numpy as np
from astropy.modeling.models import Gaussian2D
from turbustat.statistics import PowerSpectrum
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.utils.console import ProgressBar

from galaxy_params import gal_feath as gal
from paths import iram_co21_data_path

imgsize = (2056, 2056)

sigma = 10.
amp = 1.

Ngauss = 500

image = np.zeros(imgsize)

yy, xx = np.indices(imgsize)

for i in range(Ngauss):
    yc = np.random.rand() * imgsize[0]
    xc = np.random.rand() * imgsize[1]

    image += Gaussian2D(amplitude=amp, x_mean=xc, y_mean=yc,
                        x_stddev=sigma, y_stddev=sigma)(yy, xx)


pspec1 = PowerSpectrum(fits.PrimaryHDU(image))
pspec1.run(verbose=True, fit_2D=False, high_cut=10**-1.2 / u.pix,
           fit_kwargs={'brk': 0.01 / u.pix})

raw_input("?")
plt.clf()

# Get a power-spectrum that is flat on large scales. Now let's cluster the
# sources together.

image_clust = np.zeros(imgsize)

yy, xx = np.indices(imgsize)
yy -= imgsize[0] // 2
xx -= imgsize[1] // 2

for i in ProgressBar(range(Ngauss)):
    # yc = np.random.normal(loc=imgsize[0] / 2, scale=imgsize[0] * 0.2)
    # xc = np.random.normal(loc=imgsize[1] / 2, scale=imgsize[1] * 0.2)
    yc = np.random.exponential(scale=imgsize[0] * 0.2) * \
        np.random.choice([-1, 1])
    xc = np.random.exponential(scale=imgsize[1] * 0.2) * \
        np.random.choice([-1, 1])

    image_clust += Gaussian2D(amplitude=amp, x_mean=xc, y_mean=yc,
                              x_stddev=sigma, y_stddev=sigma)(yy, xx)

pspec2 = PowerSpectrum(fits.PrimaryHDU(image_clust))
pspec2.run(verbose=True, fit_2D=False, high_cut=10**-1.2 / u.pix,
           fit_kwargs={'brk': 0.01 / u.pix})

raw_input("?")
plt.clf()


# Still fairly flat at large-scale, but the largest scale now has a
# significant increase compare to the first case

# Now allow sigma to vary a bit

image_clust_sigvar = np.zeros(imgsize)

yy, xx = np.indices(imgsize)
yy -= imgsize[0] // 2
xx -= imgsize[1] // 2

sigmas = []

posns = []

for i in range(Ngauss):
    sigma_i = np.random.normal(loc=15., scale=5.)
    sigmas.append(sigma_i)

    # yc = np.random.normal(loc=imgsize[0] / 2, scale=imgsize[0] * 0.2)
    # xc = np.random.normal(loc=imgsize[1] / 2, scale=imgsize[1] * 0.2)
    yc = np.random.exponential(scale=imgsize[0] * 0.2) * \
        np.random.choice([-1, 1])
    xc = np.random.exponential(scale=imgsize[1] * 0.2) * \
        np.random.choice([-1, 1])

    posns.append([yc, xc])

    image_clust_sigvar += Gaussian2D(amplitude=amp, x_mean=xc, y_mean=yc,
                                     x_stddev=sigma_i,
                                     y_stddev=sigma_i)(yy, xx)

pspec3 = PowerSpectrum(fits.PrimaryHDU(image_clust_sigvar))
pspec3.run(verbose=True, fit_2D=False, high_cut=10**-1.2 / u.pix,
           fit_kwargs={'brk': 0.01 / u.pix})

raw_input("?")
# Now we get an upper slope of ~ -1. This version is similar in shape to the
# CO(2-1) power spectrum.

plt.close()

# Now let's take the azimuthally-averaged CO profile and
# generate a disk from it.

# exp_model = lambda r: 12 * np.exp(-r / 2.1)
exp_model = lambda r: (1 / 2.1) * np.exp(-r / 2.1)

mom0 = fits.open(iram_co21_data_path("m33.co21_iram.mom0.fits"))[0]

radii = gal.radius(header=mom0.header).to(u.kpc).value

disk_model = exp_model(radii)

pspec_disk = PowerSpectrum(fits.PrimaryHDU(disk_model))
pspec_disk.run(verbose=True, fit_2D=False, high_cut=10**-1.2 / u.pix,
               apodize_kernel='tukey')

raw_input("?")

plt.close()

plt.loglog(pspec1.freqs, pspec1.ps1D,
           label='Width set; rand. posn')
plt.loglog(pspec2.freqs, pspec2.ps1D,
           label='Width set; clust posn')
plt.loglog(pspec3.freqs, pspec3.ps1D,
           label=r'Width $\mathcal{N}(15, 5)$; clust posn')
plt.loglog(pspec_disk.freqs, pspec_disk.ps1D, label='Exp. Disk')
plt.grid()
plt.legend(frameon=True)
plt.xlabel("Freq (1 / pix)")
plt.tight_layout()
plt.savefig("/home/eric/Dropbox/Various Plots/gaussian_sources_pspec.png")
plt.close()
