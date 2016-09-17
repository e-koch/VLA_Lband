
'''
Simultaneously average all RRL SPWs together. Maybe there's a bit of signal?

Using stacking from CASA guides:
https://casaguides.nrao.edu/index.php?title=Stacking_Multiple_Spectral_Lines_at_Same_Position
'''

import os
import json

from tasks import imhead, imsmooth, impbcor, immath

root = os.path.expanduser("~/Data_1/14B-088_HRL/")

lines = ["H152alp", "H153alp", "H158alp", "H164alp", "H166alp"]
# lines = ["H152alp", "H153alp", "H166alp"]

# Load in line file information
line_file = os.path.expanduser(
    "~/Dropbox/code_development/VLA_Lband/14B-088/Lines/spw_dict.txt")
with open(line_file, 'r') as f:
    line_dict = json.load(f)
restfreq = [line_dict[line] for line in lines]


image_names = [os.path.join(root, line,
                            "dirty_images/{}_14B-088_dirty.image".format(line))
               for line in lines]

output_stacked_image = \
    os.path.join(root, "combined_imaging/HRL_14B-088_combined.image")

# Getting the sizes of the beam.
beammajor = []
beamminor = []
beampa = []

for image in image_names:
    beammajor.append(imhead(imagename=image, mode='get', hdkey='beammajor'))
    beamminor.append(imhead(imagename=image, mode='get', hdkey='beamminor'))
    beampa.append(imhead(imagename=image, mode='get', hdkey='beampa'))

# determining what size the spectral line images should be smoothed to.
targetbeammajor = max(beammajor)
targetidx = beammajor.index(targetbeammajor)
targetbeamminor = beamminor[targetidx]
targetbeampa = beampa[targetidx]

# Apply PB correction
for image in image_names:
    outputimage = image + ".pbcor"
    impbcor(imagename=image, pbimage=image.replace(".image", ".flux"),
            outfile=outputimage)

image_names = [name + ".pbcor" for name in image_names]

# smoothing the images
for image in image_names:
    outputimage = image + '.smooth'
    if image != image_names[targetidx]:
        imsmooth(imagename=image, outfile=outputimage, kernel='gauss',
                 major=str(targetbeammajor['value']) + targetbeammajor['unit'],
                 minor=str(targetbeamminor['value']) + targetbeamminor['unit'],
                 pa=str(targetbeampa['value']) + targetbeampa['unit'],
                 targetres=True)
    else:
        # need to skip the one with the largest resolution or imsmooth fails.
        os.system('cp -ir ' + image + ' ' + outputimage)

# creating the expression for immath. It doesn't appear to know about mean.
myim = ['IM' + str(i) for i in range(len(image_names))]
myexp = '(' + '+'.join(myim) + ')/' + str(float(len(image_names)))

image_names = [name + ".smooth" for name in image_names]

# combining the images
immath(imagename=image_names,
       outfile=output_stacked_image, mode='evalexpr', expr=myexp)
