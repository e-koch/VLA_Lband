
# Give a directory with a set of MSs, concatenate them, run uvcontsub,
# then image. This isn't meant to be completely automated. Parameters still
# need to be set in this script.

import os
import sys
import json

# from casa_tools import myconcat, myclean, myuvcontsub
from tasks import concat, clean, uvcontsub
from tasks import rmtables


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]


line_direc = sys.argv[-1]

line_name = line_direc.rstrip("/").split("/")[-1]
project_name = "14B-088"

# Read in rest frequency from reference saved in this repo
line_file = \
    os.path.expanduser(
        "~/Dropbox/code_development/VLA_Lband/14B-088/Lines/spw_dict.txt")
with open(line_file, 'r') as f:
    line_dict = json.load(f)
rest_freq = str(line_dict[line_name])

# Toggle on/off different operations
do_concat = False
dirty_cube_nosub = False
contsub = False
dirty_cube = True
cont_cube = False
clean_cube = False

# Parameters needed for multiple parts
concat_ms = os.path.join(line_direc, line_name + "_" + project_name + ".ms")
contsub_ms = concat_ms + ".contsub"
cont_ms = concat_ms + ".cont"
dirty_image_direc = os.path.join(line_direc, "dirty_images")
clean_image_direc = os.path.join(line_direc, "clean_images")

delete_old = True

# Contsub parameters
fitspw = "0:50~100;400~450"
excludechans = False

# General clean parameters
imsize = [2560, 2560]
cell = '3arcsec'
mode = 'velocity'
nchan = 28
width = "10.0km/s"
start = "-320.0km/s"
threshold = "1.0mJy"
field = "M33*"
# field = ",".join(["M33_{}".format(i) for i in range(1, 15) if i not in [3, 7]]) + ", M33_7_center"
phasecenter = 'J2000 01h33m50.904 +30d39m35.79'
spw = ""
imagermode = "mosaic"
multiscale = []
outframe = "LSRK"
veltype = "radio"
minpb = 0.1
weighting = 'natural'
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

    concat(vis=ms_names, concatvis=concat_ms)

if dirty_cube_nosub:

    # Check that the concatenated MS exists
    if not os.path.exists(concat_ms):
        raise IOError("Concatenated MS does not exist in the given directory.")

    if not os.path.exists(dirty_image_direc):
        os.mkdir(dirty_image_direc)

    out_image = os.path.join(dirty_image_direc,
                             line_name + "_" + project_name + "_dirty_nocontsub")

    rmtables(out_image + ".*")

    clean(vis=concat_ms, imagename=out_image, niter=0, imsize=imsize,
          cell=cell, mode=mode, nchan=nchan, width=width, start=start,
          threshold=threshold, field=field, phasecenter=phasecenter, spw=spw,
          imagermode=imagermode, multiscale=multiscale, outframe=outframe,
          veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
          restfreq=restfreq, usescratch=usescratch,
          interpolation=interpolation, pbcor=True)

if contsub:
    # Check that the concatenated MS exists
    if not os.path.exists(concat_ms):
        raise IOError("Concatenated MS does not exist in the given directory.")

    if os.path.exists(contsub_ms) and delete_old:
        rmtables(contsub_ms)

    uvcontsub(vis=concat_ms, fitspw=fitspw, excludechans=excludechans,
              fitorder=0, want_cont=False)

if dirty_cube:
    # Check that the contsub MS exists
    if not os.path.exists(contsub_ms):
        raise IOError("Contsub MS does not exist in the given directory.")

    if not os.path.exists(dirty_image_direc):
        os.mkdir(dirty_image_direc)

    out_image = os.path.join(dirty_image_direc,
                             line_name + "_" + project_name + "_dirty")
                             # line_name + "_" + project_name + "_dirty_noM33_3_equaltime")

    rmtables(out_image + ".*")

    clean(vis=contsub_ms, imagename=out_image, niter=0, imsize=imsize,
          cell=cell, mode=mode, nchan=nchan, width=width, start=start,
          threshold=threshold, field=field, phasecenter=phasecenter, spw=spw,
          imagermode=imagermode, multiscale=multiscale, outframe=outframe,
          veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
          restfreq=restfreq, usescratch=usescratch,
          interpolation=interpolation, pbcor=True, scan="0~783")

if cont_cube:
    # Check that the contsub MS exists
    if not os.path.exists(cont_ms):
        raise IOError("Continuum MS does not exist in the given directory.")

    if not os.path.exists(dirty_image_direc):
        os.mkdir(dirty_image_direc)

    out_image = os.path.join(dirty_image_direc,
                             line_name + "_" + project_name + "_contonly_dirty")

    rmtables(out_image + ".*")

    clean(vis=cont_ms, imagename=out_image, niter=0, imsize=imsize,
          cell=cell, mode=mode, nchan=nchan, width=width, start=start,
          threshold=threshold, field=field, phasecenter=phasecenter, spw=spw,
          imagermode=imagermode, multiscale=multiscale, outframe=outframe,
          veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
          restfreq=restfreq, usescratch=usescratch,
          interpolation=interpolation, pbcor=True)

if clean_cube:
    # Check that the contsub MS exists
    if not os.path.exists(concat_ms):
        raise IOError("Concatenated MS does not exist in the given directory.")

    if not os.path.exists(clean_image_direc):
        os.mkdir(clean_image_direc)

    out_image = os.path.join(clean_image_direc,
                             line_name + "_" + project_name)

    rmtables(out_image + ".*")

    clean(vis=contsub_ms, imagename=out_image, niter=10000, interactive=True,
          pbcor=True, imsize=imsize,
          cell=cell, mode=mode, nchan=nchan, width=width, start=start,
          threshold=threshold, field=field, phasecenter=phasecenter, spw=spw,
          imagermode=imagermode, multiscale=multiscale, outframe=outframe,
          veltype=veltype, minpb=0.3, weighting=weighting, robust=robust,
          restfreq=restfreq, usescratch=usescratch,
          interpolation=interpolation)
