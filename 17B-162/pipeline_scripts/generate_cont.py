'''
Create the cont.dat file for the spectral lines.

Assumed to be run on the spectral line only. So the 0th SPW is HI.
'''

import astropy.units as u
import numpy as np

vsys = -180 * u.km / u.s

bw_end_frac = 0.1

m33_mw_avoid_hi = [70, -320] * u.km / u.s
m33_avoid = [40, -40, -80, -280] * u.km / u.s

# Pol cals split out already
fields = ["0137+331=3C48", "M33_2", "M33_14", "M33_6", "M33_7_center",
          "M33_12", "M33_11", "M33_8"]

spw_dict = {0: {"obs_freq": 1419.121832 * u.MHz,
                "bwidth": 4. * u.MHz,
                "rest_freq": 1420.40580 * u.MHz},
            1: {"obs_freq": 1424.448260 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1424.73359 * u.MHz},
            2: {"obs_freq": 1477.084783 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1477.33457 * u.MHz},
            3: {"obs_freq": 1612.043787 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1612.23100 * u.MHz},
            4: {"obs_freq": 1651.372604 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1651.54111 * u.MHz},
            5: {"obs_freq": 1665.241392 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1665.40180 * u.MHz},
            6: {"obs_freq": 1667.199579 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1667.35900 * u.MHz},
            7: {"obs_freq": 1720.397384 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1720.52990 * u.MHz},
            8: {"obs_freq": 1818.166649 * u.MHz,
                "bwidth": 2. * u.MHz,
                "rest_freq": 1818.24591 * u.MHz},
            9: {"obs_freq": 1854.184798 * u.MHz,
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
