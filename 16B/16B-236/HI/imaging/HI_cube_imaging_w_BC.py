
'''
Imaging tests for combining A+B+C configs
'''

from tasks import tclean

# First do a simple 1-channel test using the four nearest B+C pointings
# centered on the A pointing

# execfile(os.path.expanduser("~/code/VLA_Lband/16B/spw_setup.py"))
execfile(os.path.expanduser("~/Dropbox/code_development/VLA_Lband/16B/spw_setup.py"))

ms_list = ["16B-236_HI_spw_0_LSRK.ms.contsub",
           "14B-088_HI_LSRK.16B236_fields.ms.contsub",
           "17B-162_HI_spw_0_LSRK.16B236_fields.ms.contsub"]

spw_num = 0

tclean(vis=ms_list,
       datacolumn='corrected',
       imagename='M33_16B-236_14B_17B_{0}_spw_{1}.dirty'.
            format(linespw_dict[spw_num][0], spw_num),
       spw="0",
       field='*',
       imsize=[10000, 10000],
       cell="0.3arcsec",
       specmode='cube',
       start="-140km/s",
       width="5km/s",
       nchan=1,
       startmodel=None,
       gridder='mosaic',
       weighting='natural',
       niter=0,
       threshold='3.2mJy/beam',
       phasecenter='J2000 01h33m33.191 +30d32m06.72',
       restfreq=linespw_dict[spw_num][1],
       outframe='LSRK',
       pblimit=0.7,
       usemask='pb',
       mask=None,
       deconvolver='hogbom',
       pbcor=False,
       veltype='radio',
       chanchunks=-1,
       restoration=False,
       # parallel=True,
       parallel=False,
       )
