
'''
Create test images of the central channels using the calibrator.
'''

import sys
import os

pipepath = '/lustre/aoc/observers/nm-7669/canfar_scripts/EVLA_pipeline1.3.0/'

execfile(pipepath+"EVLA_pipe_restore.py")

ms_name = ms_active

# Based on C config beam
myimsize = 2048
mycell = "3arcsec"

try:
    what_to_image = sys.argv[1]
except IndexError:
    what_to_image = str(raw_input("Image cals, sources or both? : "))

# if what_to_image != "cals" or what_to_image != "sources" or what_to_image != "both":
#     raise TypeError("what_to_image must be 'cals', 'sources', or 'both'.")

if not os.path.exists("test_images"):
    os.mkdir("test_images")

# Get the number of channels in each SPW
tb.open(ms_name + "/SPECTRAL_WINDOW")
nchans = tb.getcol("NUM_CHAN")
ref_freqs = tb.getcol('REF_FREQUENCY')
tb.close()

for i in spws:

    print("On SPW "+str(i)+" of "+str(len(spws)))

    if what_to_image == "cals" or what_to_image == "both":
        for cal_field in calibrator_field_list:
            print("Imaging calibrator in SPW "+str(i))

            try:
                imname = "test_images/calibrator_field" + str(cal_field) + \
                    "_spw" + str(i)

                default('clean')
                vis = 'calibrators.ms'
                imagename = imname
                outlierfile = ''
                field = str(cal_field)
                spw = str(i)
                selectdata = False
                mode = 'channel'
                nchan = 1
                width = 1
                start = 0
                reffreq = ''
                gridmode = ''
                niter = 0
                threshold = '0.0mJy'
                imagermode = ''
                imsize = [myimsize, myimsize]
                cell = [mycell, mycell]
                phasecenter = ''
                restfreq = ''
                stokes = 'I'
                weighting = 'uniform'
                uvtaper = False
                modelimage = ''
                restoringbeam = ['']
                pbcor = False
                minpb = 0.2
                calready = False
                allowchunk = False
                async = False
                clean()

            except Exception, e:
                logprint("Problem with " + str(cal_field) + " " + str(i) + " " + imname)
                print e

    # Now image the central channels of the source
    if what_to_image == "sources" or what_to_image == "both":

        print("Imaging source in SPW "+str(i))

        try:
            imname = "test_images/source_field" + "_spw" + str(i)

            default('clean')
            vis = ms_name
            imagename = imname
            outlierfile = ''
            field = "M33*"
            spw = str(i)
            selectdata = False
            mode = 'mfs'
            nterms = 1
            reffreq = ''
            gridmode = ''
            niter = 0
            threshold = '0.0mJy'
            imagermode = 'mosaic'
            imsize = [myimsize, myimsize]
            cell = [mycell, mycell]
            phasecenter = 'J2000 01h33m50.904 +30d39m35.79'
            restfreq = str(ref_freqs[i])+"Hz"
            stokes = 'I'
            weighting = 'natural'
            uvtaper = False
            pbcor = False
            minpb = 0.2
            calready = False
            allowchunk = False
            async = False
            clean()

        except Exception, e:
            logprint("Problem with " + str(i) + " " + imname)
            print e
