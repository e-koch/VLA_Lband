
#  Line SPW setup for 17B-162 w/ rest frequencies
# SPW $: [Name, Restfreq, Num chans, Chans for uvcontsub]

# Notes for uvcontsub channels:
# - Avoid first and last 5% of channels in a SPW
# - Avoid -80 to -280 km/s for M33
# - Avoid -30 to 30 km/s for MW

# Channel offset b/w cubes made with test_line_imaging.py (excludes edge
# channels) and the full SPW (since we split from the full SPW):
# Halp - 7 channels
# OH - 13 channels
# HI - 205 channels

linespw_dict = {0: ["HI", "1.420405752GHz", 4096, "1240~1560;2820~3410"],
                1: ["H166alp", "1.42473GHz", 128, "32~47;107~114"],
                2: ["H164alp", "1.47734GHz", 128, "32~47;107~114"],
                3: ["OH1612", "1.612231GHz", 256, "53~88;223~240"],
                4: ["H158alp", "1.65154GHz", 128, "32~47;107~114"],
                5: ["OH1665", "1.6654018GHz", 256, "53~88;223~240"],
                6: ["OH1667", "1.667359GHz", 256, "53~88;223~240"],
                7: ["OH1720", "1.72053GHz", 256, "53~88;223~240"],
                8: ["H153alp", "1.81825GHz", 128, "32~47;107~114"],
                9: ["H152alp", "1.85425GHz", 128, "32~47;107~114"]}
