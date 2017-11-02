
'''
Paper 2 analysis script.

Comparison between HI and CO.
'''

from os.path import join as osjoin

from paths import c_hi_analysispath, ancillary_path, root

# Masking and moment-making scripts

# Make velocity corrected cubes
execfile(osjoin(ancillary_path, "IRAM30M_CO21/14B-088/co_cube_shift.py"))

# Analysis & Figures

# Make a column density comparison plot
execfile(c_hi_analysispath("co_comparisons/co_hi_moment0_figure.py"))

# Compare velocity surfaces
execfile(c_hi_analysispath("co_comparisons/co_hi_veloffset.py"))

# Surface density profiles
execfile(c_hi_analysispath("co_comparisons/co_radial_profile.py"))

# Create stacked profiles
execfile(c_hi_analysispath("co_comparisons/total_stacked_profiles.py"))
# Analyze those stacked profiles
execfile(c_hi_analysispath("co_comparisons/stacking_analysis.py"))


# Limited spectral decomposition
execfile(c_hi_analysispath("co_comparisons/co_hi_linewidth_ratio.py"))
# Analysis of those fits
execfile(c_hi_analysispath("co_comparisons/h2_hi_ratio_analysis.py"))


# Create an HI cube matching the CO spectral resolution
execfile(osjoin(root, "14B-088/HI/imaging/HI_spectral_downsample_to_CO21.py"))
# Adaptive threshold comparison
execfile(c_hi_analysispath("co_comparisons/co_vs_hi_boundaries.py"))

# Create tables of column density
execfile(c_hi_analysispath("co_comparisons/h2_hi_ratios.py"))

# Comparisons to the Krumholz models and other population studies
execfile(c_hi_analysispath("co_comparisons/krumholz_models.py"))

# HI properties w/ and w/o CO detected
execfile(c_hi_analysispath("co_comparisons/hi_props_vs_co_detection.py"))
