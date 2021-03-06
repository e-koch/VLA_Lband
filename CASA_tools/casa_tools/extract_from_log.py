
'''
Extract info from a call to the clean task from CASA logs.
'''

import re
from itertools import izip
from datetime import datetime
from astropy import units as u
import numpy as np
from astropy.table import Table

# Define some strings for re
all_time_date = r"^[0-9]{4}-[0-9]{2}-[0-9]{2}\s[0-9]{2}:[0-9]{2}:[0-9]{2}\s"
casa_datetime_format = r'%Y-%m-%d %H:%M:%S'

info = r"INFO\s"
warn = r"WARN\s"
err = r"ERROR\s"

numbers = r"[-+]?\d*\.\d+|\d+"


def collect_clean_results(log_files, filename=None, format='ascii.csv',
                          show_in_browser=False):
    '''
    Loop through the list of given log files, extract results from the clean
    calls, and save as a csv file.

    Parameters
    ----------
    log_files : list or np.ndarray
        List or array of the log file names.
    filename : str, optional
        Name of file to save with clean results. If None is given, no file is
        saved.
    format : str of filetype
        Filetype to save the table as. See the list of writers available for
        `~astropy.table` here:
        `http://docs.astropy.org/en/stable/io/unified.html#built-in-readers-writers`_
    show_in_browser : bool, optional
        Displays the table in a web browser.

    '''
    results_dict = {"Name": [],
                    "Reached Threshold": [],
                    "Max Residual": [],
                    "Iterations": [],
                    "Time Elapsed": []}

    for i, log in enumerate(log_files):
        results = CleanResults(log)
        try:
            results.run_all()
            # Extract units
            bright_unit = results.max_residuals.unit
            time_unit = results.time_elapsed.unit

            results_dict["Name"].append(log.rstrip(".log"))
            results_dict["Reached Threshold"].append(results.finished)
            results_dict["Max Residual"].append(results.max_residuals.value)
            results_dict["Iterations"].append(results.niters)
            results_dict["Time Elapsed"].append(results.time_elapsed.value)

        except Warning as e:
            print("Failed for log: " + log)
            print(e)

            results_dict["Name"].append(log.rstrip(".log"))
            results_dict["Reached Threshold"].append(False)
            results_dict["Max Residual"].append(np.NaN)
            results_dict["Iterations"].append(0)
            results_dict["Time Elapsed"].append(np.NaN)

    # Add units back on
    results_dict["Max Residual"] *= bright_unit
    results_dict["Time Elapsed"] *= time_unit

    # Now gather into a table.
    t = Table(results_dict.values(), names=results_dict.keys())

    if filename is not None:
        t.write(filename, format=format)

    if show_in_browser:
        t.show_in_browser()


class CleanResults(object):
    """
    Read the results of running clean from a log file.

    Parameters
    ----------
    filename : str
        Name of the log file to search.
    """

    def __init__(self, filename):
        self.filename = filename

        self._lines = load_log(filename)

        self._line_ranges = None

        self._max_residuals = None

        self._time_elapsed = None

    @property
    def lines(self):
        return self._lines

    def search_log(self, expression, view=None, return_linenum=True):
        '''
        Search through the log for a given expression.
        Return the matched lines.
        '''
        re_express = re.compile(expression)

        if view is None:
            view = slice(None)

        view = fill_in_slice(view, len(self.lines))

        linenum_gen = xrange(view.start, view.stop, view.step)

        matched_lines = []
        matched_line_nums = []
        for i, line in izip(linenum_gen, self.lines[view]):
            search = re_express.search(line)
            if search:
                matched_lines.append(line)
                if return_linenum:
                    matched_line_nums.append(i)

        if not matched_lines:
            Warning("No matches found.")

        if len(matched_lines) == 1:
            if return_linenum:
                return matched_lines[0], matched_line_nums[0]
            return matched_lines[0]

        if return_linenum:
            return zip(matched_lines, matched_line_nums)
        return matched_lines

    @property
    def finished(self):
        '''
        Did CLEAN reach the given threshold?
        '''
        return self._finished_calls

    def get_finished(self):

        finish_re = all_time_date + info + \
            "*MFMSCleanImageSkyModel::solve\s*Reached*"

        if not self.line_ranges:
            self.get_line_ranges()

        if isinstance(self.line_ranges[0], int):

            start, stop = self.line_ranges
            finish_match = \
                self.search_log(finish_re, view=slice(start, stop))

            self._finished_calls = False if not finish_match else True

        else:
            finished_calls = []
            for clean_range in self.line_ranges:
                start, stop = clean_range
                finish_match = self.search_log(finish_re,
                                               view=slice(start, stop))

                if not finish_match:
                    finished_calls.append(False)
                else:
                    finished_calls.append(True)

            self._finished_calls = finished_calls

    @property
    def line_ranges(self):
        return self._line_ranges

    def get_line_ranges(self):
        '''
        Find the beginning and end of CLEAN.
        '''

        start_re = all_time_date + info + \
            "*clean::::.\s####.*Begin Task: clean.*"
        stop_re = all_time_date + info + "*clean::::\s####.*End Task: clean.*"

        start_search = self.search_log(start_re)
        if start_search:
            start_lines = start_search[1]
            self._error = False
        else:
            raise Warning("Could not find CASA clean call in log.")
            self._error = True

        stop_search = self.search_log(stop_re)
        if stop_search:
            stop_lines = stop_search[1]
            self._error = False
        else:
            Warning("Could not find end to clean call. "
                    "An error likely occurred in CASA. "
                    "Setting the end to the final log line.")
            stop_lines = len(self.lines) - 1
            self._error = True

        # If they aren't equal, there was an error (no end line)
        # Must be the last clean call, since casa always crashes
        # in my experience.
        try:
            if len(start_lines) != len(stop_lines):
                Warning("One of the CLEAN class failed.")
                self._error = True
                start_lines.pop(-1)
            self._line_ranges = zip(start_lines, stop_lines)
        except TypeError:
            self._line_ranges = [start_lines, stop_lines]
            self._error = False

    @property
    def error(self):
        return self._error

    @property
    def time_elapsed(self):
        return self._time_elapsed

    def get_time_elapsed(self, output_unit=u.min):
        '''
        Find the time needed for CLEAN to run.
        '''

        if not self.line_ranges:
            self.get_line_ranges()

        if isinstance(self.line_ranges[0], int):

            start, stop = self.line_ranges

            start_time = datetime.strptime(casa_time(self.lines[start]),
                                           casa_datetime_format)
            stop_time = datetime.strptime(casa_time(self.lines[stop]),
                                          casa_datetime_format)

            self._time_elapsed = \
                time_difference(start_time, stop_time, output_unit=output_unit)

        else:
            self._time_elapsed = []
            for clean_range in self.line_ranges:
                start, stop = clean_range

                start_time = datetime.strptime(casa_time(self.lines[start]),
                                               casa_datetime_format)
                stop_time = datetime.strptime(casa_time(self.lines[stop]),
                                              casa_datetime_format)

                diff_time = \
                    time_difference(start_time, stop_time,
                                    output_unit=output_unit)
                self._time_elapsed.append(diff_time)

    @property
    def max_residuals(self):
        return self._max_residuals

    def get_max_residuals(self):

        res_re = all_time_date + info + \
            "*MFMSCleanImageSkyModel::solve\s*Final maximum*"

        if not self.line_ranges:
            self.get_line_ranges()

        if isinstance(self.line_ranges[0], int):

            start, stop = self.line_ranges
            res_match = \
                self.search_log(res_re, view=slice(start, stop))

            if not res_match:
                Warning("Could not find final residual value.")
                self._max_residuals = np.NaN * u.Jy / u.beam
            else:
                self._max_residuals = \
                    float(re.findall(numbers, res_match[
                          0])[-1]) * u.Jy / u.beam

        else:
            self._max_residuals = []
            for clean_range in self.line_ranges:
                start, stop = clean_range
                res_match = \
                    self.search_log(res_re, view=slice(start, stop))

                if not res_match:
                    Warning("Could not find final residual value.")
                    self._max_residuals.append(np.NaN * u.Jy / u.beam)
                else:
                    residual = \
                        float(re.findall(numbers, res_match)
                              [-1]) * u.Jy / u.beam
                    self._max_residuals.append(residual)

    @property
    def niters(self):
        return self._niters

    def get_niters(self):

        iter_re = all_time_date + info + \
            "*MFMSCleanImageSkyModel::solve\s*Clean used*"

        if not self.line_ranges:
            self.get_line_ranges()

        if isinstance(self.line_ranges[0], int):

            start, stop = self.line_ranges
            iter_match = \
                self.search_log(iter_re, view=slice(start, stop),
                                return_linenum=False)

            if not iter_match:
                Warning("Could not find number of iterations used.")
                self._niters = np.NaN
            else:
                # Take the last one, since it is printed out for each
                # major cycle.
                if isinstance(iter_match, list):
                    last_match = iter_match[-1]
                else:
                    last_match = iter_match
                self._niters = \
                    int(re.findall(numbers, last_match)[-1])

        else:
            self._niters = []
            for clean_range in self.line_ranges:
                start, stop = clean_range
                iter_match = \
                    self.search_log(iter_re, view=slice(start, stop),
                                    return_linenum=False)

                if not iter_match:
                    Warning("Could not find number of iterations used.")
                    self._niters.append(np.NaN)
                else:
                    if isinstance(iter_match, list):
                        last_match = iter_match[-1]
                    else:
                        last_match = iter_match

                    iters = \
                        int(re.findall(numbers, last_match)[-1])
                    self._max_residuals.append(iters)

    def run_all(self, time_unit=u.min):

        self.get_line_ranges()
        self.get_finished()
        self.get_max_residuals()
        self.get_time_elapsed(output_unit=time_unit)
        self.get_niters()

    def info_dict(self):
        if isinstance(self.line_ranges[0], int):
            return {"Finished": self.finished,
                    "Max Residual": self.max_residuals,
                    "Time Elapsed": self.time_elapsed,
                    "Iterations": self.niters}
        else:
            results_dicts = []
            for i in xrange(len(self.line_ranges[0])):
                results_dicts.append(
                    {"Finished": self.finished[i],
                     "Max Residual": self.max_residuals[i],
                     "Time Elapsed": self.time_elapsed[i],
                     "Iterations": self.niters[i]})
            return results_dicts

    def __repr__(self):
        if isinstance(self.line_ranges[0], int):
            return "Finished: " + str(self.finished) + "\nMax Residual: " + \
                str(self.max_residuals) + "\nIterations: " + \
                str(self.niters) + "\nTime Elapsed: " + \
                str(self.time_elapsed.round(3))
        else:
            for i in xrange(len(self.line_ranges[0])):
                return "Clean " + str(i + 1) + " Finished: " + \
                    str(self.finished[i]) + "\n  Max Residual: " + \
                    str(self.max_residuals[i]) + "\n  Iterations: " + \
                    str(self.niters[i]) + "\n  Time Elapsed: " + \
                    str(self.time_elapsed[i].round(3))


def load_log(logfile):
    '''
    Load the lines of a log file in.
    '''

    with open(logfile) as f:
        lines = f.readlines()

    return lines


def fill_in_slice(view, list_len):
    '''
    To keep track of lines in the log, fill in
    undefined slice parameters with defaults.
    '''

    if not view.start:
        start = 0
    else:
        start = view.start

    if not view.stop:
        stop = list_len
    else:
        stop = view.stop

    if not view.step:
        step = 1
    else:
        step = view.step

    return slice(start, stop, step)


def casa_time(line):
    return line.split("\t")[0]


def time_difference(time1, time2, output_unit=u.min):

    diff = time2 - time1

    seconds_diff = diff.total_seconds() * u.s

    return seconds_diff.to(output_unit)
