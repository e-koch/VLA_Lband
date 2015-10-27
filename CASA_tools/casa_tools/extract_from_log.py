
'''
Extract info from the CASA logs.
'''

import re
from itertools import izip

# Define some strings for re
all_time_date = "^[0-9]{4}-[0-9]{2}-[0-9]{2}\s[0-9]{2}:[0-9]{2}:[0-9]{2}\s"

info = "INFO\s"
warn = "WARN\s"
err = "ERROR\s"


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
                return matched_lines[0], matched_line_nums[1]
            return matched_lines[0]

        if return_linenum:
            return zip(matched_lines, matched_line_nums)
        return matched_lines

    @property
    def finished(self):
        '''
        Did CLEAN reach the given threshold?
        '''
        return self._finished

    def _finished(self, view=None):

        finish_re = all_time_date+info+"*MFMSCleanImageSkyModel::solve\s*Reached*"

        finish_match = self.search_log(finish_re, view=view)


        if isinstance(finish_match, list):
            return [True] * len(finish_match)
        else:
            pass

    @property
    def line_ranges(self):
        return self._line_ranges

    def start_stop_range(self):
        '''
        Find the beginning and end of CLEAN.
        '''

        start_re = all_time_date+info+"*clean::::\s####.*End Task: clean.*"
        stop_re = all_time_date+info+"*clean::::.\s####.*Begin Task: clean.*"

        start_lines = self.search_log(start_re)[1]
        stop_lines = self.search_log(stop_re)[1]

        # If they aren't equal, there was an error (no end line)
        # Must be the last clean call, since casa always crashes
        # in my experience.
        if len(start_lines) != len(stop_lines):
            Warning("One of the CLEAN class failed.")
            self._error = True
            start_lines.pop(-1)

        return zip(start_lines, stop_lines)

    @property
    def error(self):
        return self._error

    def max_residual(self):
        pass

    def time_elapsed():
        pass

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
