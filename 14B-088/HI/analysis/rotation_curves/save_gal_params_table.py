
'''
Save the galactic parameters table as a latex table in the paper repo
'''

import os
import sys
from pandas import read_csv

from analysis.paths import (fourteenB_HI_data_path, paper1_tables_path)

try:
    folder_name = sys.argv[1]
except IndexError:

    diskfit_runs = [dir for dir in os.listdir(fourteenB_HI_data_path(""))
                    if dir.startswith("diskfit")]

    print("Available diskfit runs: " + str(diskfit_runs))

    folder_name = raw_input("Give folder name of the diskfit run: ")

params = os.path.basename(folder_name.rstrip("/"))

param_name = \
    fourteenB_HI_data_path("{}/rad.out.params.csv".format(folder_name))

param_table = read_csv(param_name)

param_table.to_latex(paper1_tables_path("galactic_parameters_{}.tex".format(params)))
