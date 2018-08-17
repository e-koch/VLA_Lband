'''
Create the cont.dat file for the spectral lines.

Assumed to be run on the spectral line only. So the 0th SPW is HI.
'''

import astropy.units as u
import numpy as np

vsys = -120 * u.km / u.s

bw_end_frac = 0.1

m33_mw_avoid_hi = [40, -60, -80, -160] * u.km / u.s
m33_avoid = [40, -40, -80, -280] * u.km / u.s

# Pol cals split out already
fields = ["3C48", "M33_Sarm"]

spw_dict = {0: {"obs_freq": 1419.3875 * u.MHz,
                "bwidth": 4. * u.MHz,
                "rest_freq": 1420.40580 * u.MHz},
            1: {"obs_freq": 1424.715 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1424.73359 * u.MHz},
            2: {"obs_freq": 1477.362 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1477.33457 * u.MHz},
            3: {"obs_freq": 1612.346 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1612.23100 * u.MHz},
            4: {"obs_freq": 1651.682 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1651.54111 * u.MHz},
            5: {"obs_freq": 1665.553 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1665.40180 * u.MHz},
            6: {"obs_freq": 1667.512 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1667.35900 * u.MHz},
            7: {"obs_freq": 1720.720 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1720.52990 * u.MHz},
            8: {"obs_freq": 1818.507 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1818.24591 * u.MHz},
            9: {"obs_freq": 1854.532 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1854.25027 * u.MHz}}


def generate_cont_regions(info_dict, vels):

    obs_freq = info_dict['obs_freq']
    bwidth = info_dict['bwidth']
    rest_freq = info_dict['rest_freq']

    min_freq = obs_freq + bw_end_frac * bwidth
    max_freq = obs_freq + (1.0 - bw_end_frac) * bwidth

    print("Start: {}".format(np.round(min_freq.to(u.GHz), 4)))

    for vel in vels:

        freq_i = vel.to(u.MHz, u.doppler_radio(rest_freq))

        print(vel, np.round(freq_i.to(u.GHz), 4))

    print("End: {}".format(np.round(max_freq.to(u.GHz), 4)))


print("Spectral Window: 0")
generate_cont_regions(spw_dict[0], m33_mw_avoid_hi)

for i in range(1, 10):
    print("Spectral Window: {}".format(i))
    generate_cont_regions(spw_dict[i], m33_avoid)
