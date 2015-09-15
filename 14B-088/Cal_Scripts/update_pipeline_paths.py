
'''
Update EVLA pipeline variables to the current system.
'''


def update_paths(pipe_dict, ms_path=None, pipepath=None):

    if ms_path is not None:
        pipe_dict['ms_active'] = ms_path
        pipe_dict['SDM_name'] = ms_path[:-3] # Cutoff '.ms'

    if pipepath is not None:
        pipe_dict['pipepath'] = pipepath

    return pipe_dict


if __name__ == '__main__':

    import sys

    pipe_var_file = str(sys.argv[-3])

    ms_path = str(sys.argv[-2])
    if ms_path == ".":
        ms_path = None

    pipepath = str(sys.argv[-1])
    if pipepath == ".":
        pipepath = None

    import shelve

    pipe_dict = shelve.open(pipe_var_file, writeback=True)

    pipe_dict = update_paths(pipe_dict, ms_path=ms_path, pipepath=pipepath)

    pipe_dict.sync()

    print "Checking!"
    print pipe_dict['ms_active']
    print pipe_dict['SDM_name']
    print pipe_dict['pipepath']

    pipe_dict.close()
