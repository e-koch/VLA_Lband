
# Give a directory with a set of MSs, concatenate them, run uvcontsub,
# then image. This isn't meant to be completely automated. Parameters still
# need to be set in this script.

import os
import sys
import json

from casa_tools import myconcat, myclean, myuvcontsub
from tasks import rmtables


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


line_direc = sys.argv[-1]

line_name = line_direc.rstrip("/").split("/")[-1]
project_name = "14B-088"

# Read in rest frequency from reference saved in this repo
line_file = \
    os.path.join(os.path.expanduser("~"),
                 "Dropbox/code_development/VLA_Lband/14B-088/Lines/spw_dict.txt")
with open(line_file, 'r') as f:
    line_dict = json.load(f)
rest_freq = str(line_dict[line_name])

# Toggle on/off different operations
do_concat = False
dirty_cube_nosub = False
contsub = True
dirty_cube = True
clean_cube = True

# Parameters needed for multiple parts
concat_ms = os.path.join(line_direc, line_name+"_"+project_name+".ms")
contsub_ms = concat_ms+".contsub"
dirty_image_direc = os.path.join(line_direc, "dirty_images")
clean_image_direc = os.path.join(line_direc, "clean_images")

delete_old = True

# Contsub parameters
fitspw = "0:50~100;400~450"
excludechans = False

# General clean parameters
imsize = [500, 500]  # [2560, 2560]
cell = '3arcsec'
mode = 'velocity'
nchan = 11
width = "20.0km/s"
start = "-290.0km/s"
thresh = "1.0mJy"
field = "M33_3"  # "M33*"
phasecenter = ''  # 'J2000 01h33m50.904 +30d39m35.79'
spw = ""
imagermode = "csclean"  # "mosaic"
multiscale = []
outframe = "LSRK"
veltype = "radio"
minpb = 0.3
weighting = 'briggs'
robust = 0.0
restfreq = rest_freq
usescratch = False
interpolation = 'linear'

if do_concat:
    ms_names = []
    for f in listdir_fullpath(line_direc):
        if f.endswith(".ms"):
            ms_names.append(f)

    rmtables(concat_ms)

    myconcat(vis=ms_names, concatvis=concat_ms)

if dirty_cube_nosub:

    # Check that the concatenated MS exists
    if not os.path.exists(concat_ms):
        raise IOError("Concatenated MS does not exist in the given directory.")

    if not os.path.exists(dirty_image_direc):
        os.mkdir(dirty_image_direc)

    out_image = os.path.join(dirty_image_direc,
                             line_name+"_"+project_name+"_dirty_nocontsub")

    rmtables(out_image+".*")

    myclean(vis=concat_ms, imagename=out_image, niter=0, imsize=imsize,
            cell=cell, mode=mode, nchan=nchan, width=width, start=start,
            thresh=thresh, field=field, phasecenter=phasecenter, spw=spw,
            imagermode=imagermode, multiscale=multiscale, outframe=outframe,
            veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
            restfreq=restfreq, usescratch=usescratch,
            interpolation=interpolation, pbcor=False)

if contsub:
    # Check that the concatenated MS exists
    if not os.path.exists(concat_ms):
        raise IOError("Concatenated MS does not exist in the given directory.")

    if os.path.exists(contsub_ms) and delete_old:
            rmtables(contsub_ms)

    myuvcontsub(vis=concat_ms, fitspw=fitspw, excludechans=excludechans,
                fitorder=0, want_cont=True)

if dirty_cube:
    # Check that the contsub MS exists
    if not os.path.exists(concat_ms):
        raise IOError("Concatenated MS does not exist in the given directory.")

    if not os.path.exists(dirty_image_direc):
        os.mkdir(dirty_image_direc)

    out_image = os.path.join(dirty_image_direc,
                             line_name+"_"+project_name+"_dirty")

    rmtables(out_image+".*")

    myclean(vis=contsub_ms, imagename=out_image, niter=0, imsize=imsize,
            cell=cell, mode=mode, nchan=nchan, width=width, start=start,
            thresh=thresh, field=field, phasecenter=phasecenter, spw=spw,
            imagermode=imagermode, multiscale=multiscale, outframe=outframe,
            veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
            restfreq=restfreq, usescratch=usescratch,
            interpolation=interpolation, pbcor=False)

if clean_cube:
    # Check that the contsub MS exists
    if not os.path.exists(concat_ms):
        raise IOError("Concatenated MS does not exist in the given directory.")

    if not os.path.exists(clean_image_direc):
        os.mkdir(clean_image_direc)

    out_image = os.path.join(clean_image_direc,
                             line_name+"_"+project_name)

    rmtables(out_image+".*")

    myclean(vis=contsub_ms, imagename=out_image, niter=1e5, interactive=True,
            pbcor=True, imsize=imsize,
            cell=cell, mode=mode, nchan=nchan, width=width, start=start,
            thresh=thresh, field=field, phasecenter=phasecenter, spw=spw,
            imagermode=imagermode, multiscale=multiscale, outframe=outframe,
            veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
            restfreq=restfreq, usescratch=usescratch,
            interpolation=interpolation)
