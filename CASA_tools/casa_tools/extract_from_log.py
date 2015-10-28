
'''
Extract info from the CASA logs.
'''

import re
from itertools import izip
from datetime import datetime
from astropy import units as u
import numpy as np

# Define some strings for re
all_time_date = r"^[0-9]{4}-[0-9]{2}-[0-9]{2}\s[0-9]{2}:[0-9]{2}:[0-9]{2}\s"
casa_datetime_format = r'%Y-%m-%d %H:%M:%S'

info = r"INFO\s"
warn = r"WARN\s"
err = r"ERROR\s"

numbers = r"[-+]?\d*\.\d+|\d+"


def load_log(logfile):
    '''
    Load the lines of a log file in.
    '''

    with open(logfile) as f:
        lines = f.readlines()

    return lines


class CleanResults(object):
    """
    Read the results of running clean from a log file.
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

        finish_re = all_time_date+info+"*MFMSCleanImageSkyModel::solve\s*Reached*"

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

        start_re = all_time_date+info+"*clean::::.\s####.*Begin Task: clean.*"
        stop_re = all_time_date+info+"*clean::::\s####.*End Task: clean.*"

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
            stop_lines = len(self.lines)-1
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

    def get_time_elapsed(self, output='minutes'):
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
                time_difference(start_time, stop_time, output=output)

        else:
            self._time_elapsed = []
            for clean_range in self.line_ranges:
                start, stop = clean_range

                start_time = datetime.strptime(casa_time(self.lines[start]),
                                               casa_datetime_format)
                stop_time = datetime.strptime(casa_time(self.lines[stop]),
                                              casa_datetime_format)

                diff_time = \
                    time_difference(start_time, stop_time, output=output)
                self._time_elapsed.append(diff_time)

    @property
    def max_residuals(self):
        return self._max_residuals

    def get_max_residuals(self):

        res_re = all_time_date+info+"*MFMSCleanImageSkyModel::solve\s*Final maximum*"

        if not self.line_ranges:
            self.get_line_ranges()

        if isinstance(self.line_ranges[0], int):

            start, stop = self.line_ranges
            res_match = \
                self.search_log(res_re, view=slice(start, stop))

            if not res_match:
                Warning("Could not find final residual value.")
                self._max_residuals = np.NaN
            else:
                self._max_residuals = \
                    float(re.findall(numbers, res_match[0])[-1]) * u.Jy/u.beam

        else:
            self._max_residuals = []
            for clean_range in self.line_ranges:
                start, stop = clean_range
                res_match = \
                    self.search_log(res_re, view=slice(start, stop))

                if not res_match:
                    Warning("Could not find final residual value.")
                    self._max_residuals.append(np.NaN)
                else:
                    residual = \
                        float(re.findall(numbers, res_match)[-1]) * u.Jy/u.beam
                    self._max_residuals.append(residual)

    def run_all(self, time_output="minutes"):

        self.get_line_ranges()
        self.get_finished()
        self.get_max_residuals()
        self.get_time_elapsed(output=time_output)

    def info_dict(self):
        if isinstance(self.line_ranges[0], int):
            return {"Finished": self.finished,
                    "Max Residual": self.max_residuals,
                    "Time Elapsed": self.time_elapsed}
        else:
            results_dicts = []
            for i in xrange(len(self.line_ranges[0])):
                results_dicts.append(
                    {"Finished": self.finished[i],
                     "Max Residual": self.max_residuals[i],
                     "Time Elapsed": self.time_elapsed[i]})
            return results_dicts

    def __repr__(self):
        if isinstance(self.line_ranges[0], int):
            return "Finished: "+str(self.finished)+"\nMax Residual: " + \
                  str(self.max_residuals)+"\nTime Elapsed: " + \
                  str(self.time_elapsed.round(3))
        else:
            for i in xrange(len(self.line_ranges[0])):
                return "Clean "+str(i+1)+" Finished: " + \
                    str(self.finished[i])+"\n  Max Residual: " + \
                    str(self.max_residuals[i])+"\n  Time Elapsed: " + \
                    str(self.time_elapsed[i].round(3))


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


def time_difference(time1, time2, output="seconds"):

    diff = time2 - time1

    if output == "seconds":
        return diff.total_seconds() * u.s
    elif output == "minutes":
        return diff.total_seconds()/60. * u.min
    elif output == "hours":
        return diff.total_seconds()/3600. * u.hour
    elif output == "days":
        return diff.total_seconds()/(3600.*24.) * u.day
    else:
        raise TypeError("output must be 'seconds', 'minutes',"
                        " 'hours', or 'days'.")
