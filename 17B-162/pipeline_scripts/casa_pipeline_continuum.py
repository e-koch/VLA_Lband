
import sys

mySDM = sys.argv[-1]

__rethrow_casa_exceptions = True
context = h_init()
context.set_state('ProjectSummary', 'observatory', 'Karl G. Jansky Very Large Array')
context.set_state('ProjectSummary', 'telescope', 'EVLA')
try:
    hifv_importdata(ocorr_mode='co', nocopy=False, vis=[mySDM],
                    createmms='automatic', asis='Receiver CalAtmosphere',
                    overwrite=False)
    # hifv_hanning(pipelinemode="automatic")
    hifv_flagdata(intents='*POINTING*,*FOCUS*,*ATMOSPHERE*,*SIDEBAND_RATIO*, \
                  *UNKNOWN*, *SYSTEM_CONFIGURATION*, *UNSPECIFIED#UNSPECIFIED*',
                  flagbackup=False, scan=True, baseband=True, clip=True,
                  autocorr=True,
                  hm_tbuff='1.5int', template=True, online=True, tbuff=0.0,
                  fracspw=0.05, shadow=True, quack=True, edgespw=True)
    hifv_vlasetjy(fluxdensity=-1, scalebychan=True, reffreq='1GHz', spix=0)
    hifv_priorcals(tecmaps=False)
    hifv_testBPdcals(weakbp=False)
    hifv_flagbaddef(pipelinemode="automatic")
    hifv_checkflag(pipelinemode="automatic")
    hifv_semiFinalBPdcals(weakbp=False)
    hifv_checkflag(checkflagmode='semi')
    hifv_semiFinalBPdcals(weakbp=False)
    hifv_solint(pipelinemode="automatic")
    hifv_fluxboot(pipelinemode="automatic")
    hifv_finalcals(weakbp=False)
    hifv_applycals(flagdetailedsum=True, flagbackup=True, calwt=[True],
                   flagsum=True, gainmap=False)
    hifv_targetflag(intents='*CALIBRATE*,*TARGET*')
    hifv_statwt(pipelinemode="automatic")
    hifv_plotsummary(pipelinemode="automatic")
    hif_makeimlist(nchan=-1, calmaxpix=300, intent='PHASE,BANDPASS')
    hif_makeimages(tlimit=2.0, hm_negativethreshold=-999.0,
                   subcontms=False, hm_masking='none', masklimit=4,
                   maxncleans=1, hm_growiterations=-999, cleancontranges=False,
                   noise='1.0Jy', hm_minbeamfrac=-999.0, target_list={}, robust=-999.0,
                   parallel='automatic', weighting='briggs', hm_noisethreshold=-999.0,
                   hm_lownoisethreshold=-999.0, npixels=0, hm_sidelobethreshold=-999.0)
finally:
    h_save()