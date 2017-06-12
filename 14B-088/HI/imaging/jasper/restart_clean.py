
import numpy as np
import os
import pandas as pd
import sys

# Give path to csv filemad from HI_get_clean_results.py
tab_name = sys.argv[1]

tab = pd.read_csv(tab_name)

# Find runs w/ finite max residuals (valid mask) that haven't reached the
# threshold
bad_chans = ~tab["Reached Threshold"] & np.isfinite(tab["Max Residual"])

if len(np.where(bad_chans)[0]) == 0:
    print("No failed cleans found!")

# Make assumption that they have cleaned too far and diverged.
for val in tab["Name"].ix[bad_chans]:
    channel_name = val.split("/")[-2]
    os.chdir(channel_name)
    # Remove clean products, job output scripts, and CASA logs
    os.system("rm -rf 14B-088_HI_LSRK.ms.contsub*.clean.*")
    os.system("rm *.log *.sub.* *.last")
    print("Done cleaning up {}".format(channel_name))

    # Now re-submit
    os.system("qsub {}.sub".format(channel_name))
    os.chdir("..")
    print("Submitted {}".format(channel_name))
