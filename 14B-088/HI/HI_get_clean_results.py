
'''
Extracts info from the CLEAN logs in separate channels.
'''

import sys
import glob
from astropy.table import Table

from casa_tools import CleanResults

def collect_clean_results(log_files, filename=None, format='ascii.csv',
                          show_in_browser=False):
    '''
    Loop through the list of given log files, extract results from the clean
    calls, and save as a csv file.
    '''
    results_dict = {"Name": [],
                    "Reached Threshold": [],
                    "Max Residual": [],
                    "Iterations": [],
                    "Time Elapsed": []}

    for log in log_files:
        results = CleanResults(log)
        results.run_all()

        results_dict["Name"].append(log.rstrip(".log"))
        results_dict["Reached Threshold"].append(results.finished)
        results_dict["Max Residual"].append(results.max_residuals)
        results_dict["Iterations"].append(results.niters)
        results_dict["Time Elapsed"].append(results.time_elapsed)

    # Now gather into a table.
    t = Table(results_dict.values(), names=results_dict.keys())

    if filename is not None:
        t.write(filename, format=format)

    if show_in_browser:
        t.show_in_browser()
