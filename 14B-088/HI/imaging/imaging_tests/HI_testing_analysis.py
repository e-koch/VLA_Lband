
'''
Split from HI_testing_comparison because I didn't feel like getting statsmodels
to play nice within CASA.
'''

import os
import pandas as pd
import statsmodels.formula.api as sm


data_path = os.path.expanduser("~/MyRAID/M33/14B-088/HI/imaging/testing")

data = pd.read_csv(os.path.join(data_path, "property_values.csv"))

# So my path to CASA 4.7 was wrong :) They all failed.
data = data[data["CASAVer"] != 470]

# Drop any other NaNs
data = data.dropna()

# Let's model some stuff! Treat the different CASA versions as random mixed
# effect

data["CASAVer"][data["CASAVer"] == 440] = 0
data["CASAVer"][data["CASAVer"] == 453] = 1
data["CASAVer"][data["CASAVer"] == 460] = 2
# data["CASAVer"][data["CASAVer"] == 470] = 3

# Create a version without any diverging cleans
good_data = data[data["peak_res"] < 0.01]

# Sum
sum_model = sm.mixedlm("sum ~ Tclean*AllFields*MScale*Mask*Model", data=data,
                       groups=data["CASAVer"]).fit(reml=False)
print(sum_model.summary())

# Can't use Tclean. Makes matrix singular.
sum_model_good = sm.mixedlm("sum ~ AllFields*MScale*Mask*Model", data=good_data,
                            groups=good_data["CASAVer"]).fit(reml=False)
print(sum_model_good.summary())
# Dominated by model (duh)

# Median
median_model = \
    sm.mixedlm("median ~ Tclean*AllFields*MScale*Mask*Model", data=data,
               groups=data["CASAVer"]).fit(reml=False)
print(median_model.summary())

# Can't use Tclean. Makes matrix singular.
median_model_good = \
    sm.mixedlm("median ~ AllFields*MScale*Mask*Model", data=good_data,
               groups=good_data["CASAVer"]).fit(reml=False)
print(median_model_good.summary())
# Dominated by inclusion of model (duh). Some interaction with AllFields

# Std
std_model = \
    sm.mixedlm("std ~ Tclean*AllFields*MScale*Mask*Model", data=data,
               groups=data["CASAVer"]).fit(reml=False)
print(std_model.summary())

# Can't use Tclean. Makes matrix singular.
std_model_good = \
    sm.mixedlm("std ~ AllFields*MScale*Mask*Model", data=good_data,
               groups=good_data["CASAVer"]).fit(reml=False)
print(std_model_good.summary())
# High significance of model

# Peak Residual
peakres_model = \
    sm.mixedlm("peak_res ~ Tclean*AllFields*MScale*Mask*Model", data=data,
               groups=data["CASAVer"]).fit(reml=False)
print(peakres_model.summary())

# Can't use Tclean. Makes matrix singular.
peakres_model_good = \
    sm.mixedlm("peak_res ~ AllFields*MScale*Mask*Model", data=good_data,
               groups=good_data["CASAVer"]).fit(reml=False)
print(peakres_model_good.summary())
# High significance of Mscale. Model also somewhat significant
# Interaction between the two also significant

# In all cases, the version of CASA used makes no difference. Excellent!
