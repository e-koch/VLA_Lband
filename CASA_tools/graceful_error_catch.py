
class CASAError(Exception):
    """
    Catch CASA errors to react appropriately in python scripts.
    Errors are already outputted to the terminal or log. This class
    only return the task where the error occurred.

    Example
    -------

    When CASA tasks fail, they return 'False'. Successful completion has no
    return, giving 'None'.

    Try imaging space potatoes:
    >>> clean_out = True if clean(vis='potato', cell='9carrots',
                                  phasecenter="SPACE",
                                  imsize=-np.inf) is None else False
    >>> if not clean_out:
    >>>     raise CASAError("Failed on clean")

    """
    pass


def catch_fail(task, **kwargs):

    # Try for casa tasks only
    try:
        task_params = task.parameters.keys()
    except AttributeError:
        task_params = task.func_code.co_varnames

    input_params = {key: value for key, value in kwargs.iteritems()
                    if key in task_params}

    if task(**input_params) is False:
        raise CASAError("Failed on " + task.__name__)
