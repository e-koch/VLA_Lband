######################################################################

# QA IMAGING OF CALIBRATORS

# NB: in development, do not use!

logprint("Starting EVLA_pipe_makeimages.py", logfileout='logs/makeimages.log')
time_list = runtiming('makeimages', 'start')

logprint("Making QA calibrator images, one per field/spw with long solint",
         logfileout='logs/makeimages.log')

default('applycal')
vis = 'calibrators.ms'
field = ''
spw = ''
intent = ''
selectdata = False
gaintable = ['finalampgaincal.g', 'finalphasegaincal.g']
gainfield = ['']
interp = ['']
spwmap = []
parang = False
calwt = False
flagbackup = True
async = False
applycal()


# Find max baseline length

ms.open('calibrators.ms')
uv_range = ms.range(["uvdist"])
uv_max = uv_range['uvdist'][1]
ms.close()

c = 2.997925e8


for jj in calibrator_field_list:
    for ii in field_spws[jj]:

        try:
            # set cell size to 1/(5.*Bmax)
            wave = c / reference_frequencies[ii]
            cellsize = 206265. * wave / uv_max / 5.
            mycell = str(cellsize) + 'arcsec'
            fwhm = 206265. * wave / 25.0
#            myimsize=getOptimumSize(fwhm/cellsize)
#            myimsize=int(round(myimsize/2.)*2)
            myimsize = 320
            imname = "field" + str(jj) + "_spw" + str(ii) + "_longsolint"

            default('clean')
            vis = 'calibrators.ms'
            imagename = imname
            outlierfile = ''
            field = str(jj)
            spw = str(ii)
            selectdata = False
            mode = 'mfs'
            nterms = 1
            reffreq = ''
            gridmode = ''
            niter = 0
            gain = 0.1
            threshold = '0.0mJy'
            psfmode = 'clark'
            imagermode = ''
            multiscale = []
            mask = []
            interactive = False
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
            logprint("Problem with " + str(jj) + " " + str(ii) + " " + imname)
            print e

logprint("Making QA calibrator images, one per field/spw with short solint",
         logfileout='logs/makeimages.log')


default('applycal')
vis = 'calibrators.ms'
field = ''
spw = ''
intent = ''
selectdata = False
gaintable = ['finalampgaincal.g', 'phaseshortgaincal.g']
gainfield = ['']
interp = ['']
spwmap = []
parang = False
calwt = False
flagbackup = True
async = False
applycal()

    logprint("Field " + str(jj) + " complete",
             logfileout='logs/makeimages.log')


logprint("Finished EVLA_pipe_makeimages.py", logfileout='logs/makeimages.log')
time_list = runtiming('makeimages', 'end')

pipeline_save()
