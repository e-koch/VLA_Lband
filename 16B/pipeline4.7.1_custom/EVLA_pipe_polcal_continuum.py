
'''
Perform the polarization calibration for the continuum SPWs
'''

import re
import os
import copy
import numpy as np

from tasks import applycal, gaincal, setjy, plotms, plotcal, polcal


def pol_fit(coeff_type='angle', reffreq=1.0, order=1):
    '''
    Fit a line to the Perley-Butler pol fractions and angles in L-band
    '''

    freqs = np.array([1.05, 1.45, 1.64, 1.95])
    polfracs = np.array([5.6, 7.5, 8.4, 9.0]) / 100.
    # polangs = np.array([-14., -11., -10., -10.])

    # The angles don't change much, and are basically constant above 1.5
    if coeff_type == "angle":
        if reffreq > 1.5:
            return np.array([-10, 0])
        else:
            # Slope from two points.
            # CASA wants solutions as a function of (freq - reffreq) / reffreq
            unscaled_slope = (-11 + 14) / (1.45 - 1.05)
            slope = reffreq * unscaled_slope
            # Intercept is defined at the reffreq
            intercept = -14 + unscaled_slope * (reffreq - 1.05)
            return np.array([intercept, slope])
    elif coeff_type == "frac":
        coeffs = np.polyfit((freqs - 1.05) / 1.05, polfracs, order)
        return coeffs[::-1]
    else:
        raise ValueError("coeffs must be 'angle' or 'frac'.")


# Gain cal for pol cals
gaincal(vis=ms_active, caltable='J0319+4130_gaincal.g', field='J0319+4130',
        spw='', solint='inf', refant=refAnt, gaintype='G', calmode='ap',
        gaintable=["final_caltables/finaldelay.k",
                   "final_caltables/finalBPcal.b"])

gaincal(vis=ms_active, caltable='3C138_gaincal.g', field='3C138',
        spw='', solint='inf', refant=refAnt, gaintype='G', calmode='ap',
        gaintable=["final_caltables/finaldelay.k",
                   "final_caltables/finalBPcal.b"])

# setjy for cross-hands based on 3C138
setjy_138 = setjy(vis=ms_active, field='3C138', standard='Perley-Butler 2013',
                  model='3C138_L.im', usescratch=False, scalebychan=True, spw='')

# It doesn't appear that we can capture the relevant fitting range with the
# setjy output. So, make the user input the values from the logger.
# inp_1 = raw_input("Input I values fit to by setjy (first and last): ")
# intensities = [float(i) for i in re.split(" |,", inp_1)]
intensities = [9.53007, 7.13258]

if len(intensities) != 2:
    raise ValueError("intensities should have two values."
                     " Got: {}".format(intensities))

# inp_2 = raw_input("Input frequencies (MHz) of given intensities: ")
# freqs = [float(i) for i in re.split(" |,", inp_2)]
freqs = [1243.45, 2010.45]

if len(freqs) != 2:
    raise ValueError("freqs should have two values."
                     " Got: {}".format(freqs))\

# Assume single powerlaw solution for the flux.
alpha = np.log(intensities[1] / intensities[0]) / np.log(freqs[1] / freqs[0])
i0 = intensities[0]

spw_reffreqs = {"0": "0.988GHz", "1": "1.116GHz", "2": "1.244GHz",
                "3": "1.372GHz", "4": "1.5GHz", "5": "1.628GHz",
                "6": "1.756GHz", "7": "1.884GHz"}

# Solve per SPW based on Perley-Butler 2013 values for 3C138.
for key in setjy_138['4']:
    if key == 'fieldName':
        continue

    reffreq = spw_reffreqs[key]
    reffreq_float = float(reffreq[:-3])

    i0 = setjy_138['4'][key]['fluxd'][0]

    d_coeff = np.deg2rad(pol_fit(coeff_type='angle', reffreq=reffreq_float))
    c_coeff = pol_fit(coeff_type='frac', reffreq=reffreq_float, order=2)

    setjy(vis=ms_active, field='3C138', standard='manual',
          spw=key, fluxdensity=[i0, 0, 0, 0], spix=[alpha, 0],
          reffreq=reffreq,
          polindex=c_coeff, polangle=d_coeff,
          scalebychan=True, usescratch=False)

default('plotms')

plotms(vis=ms_active, field='3C138', correlation='RR',
       antenna='ea01&ea02',
       xaxis='channel', yaxis='amp', ydatacolumn='model',
       plotfile='plotms_3c138-model-amp-RR.png', overwrite=True)

plotms(vis=ms_active, field='3C138', correlation='RL',
       antenna='ea01&ea02',
       xaxis='channel', yaxis='amp', ydatacolumn='model',
       plotfile='plotms_3c138-model-amp-RL.png', overwrite=True)

plotms(vis=ms_active, field='3C138', correlation='RL',
       antenna='ea01&ea02',
       xaxis='channel', yaxis='phase', ydatacolumn='model',
       plotrange=[-1, -1, -180, 180],
       plotfile='plotms_3c138-model-phase-RL.png',
       overwrite=True)

# Solve for cross-hand delays
gaincal(vis=ms_active, caltable='3C138.Kcross',
        field='3C138', spw='',
        gaintype='KCROSS', solint='inf', combine='scan', refant=refAnt,
        gaintable=["final_caltables/finaldelay.k",
                   "final_caltables/finalBPcal.b",
                   "3C138_gaincal.g"],
        gainfield=[],
        parang=True)

plotcal(caltable='3C138.Kcross', xaxis='antenna', yaxis='delay',
        figfile='plotcal_3c138-Kcross-delay.png')

# Solve for leakage terms
polcal(vis=ms_active, caltable='J0319+4130.D1',
       field='J0319+4130', spw='',
       refant=refAnt, poltype='Df', solint='inf', combine='scan',
       gaintable=["final_caltables/finaldelay.k",
                  "final_caltables/finalBPcal.b",
                  "J0319+4130_gaincal.g",
                  '3C138.Kcross'],
       gainfield=[])

plotcal(caltable='J0319+4130.D1', xaxis='chan', yaxis='amp',
        spw='', field='', iteration='antenna')
#
plotcal(caltable='J0319+4130.D1', xaxis='chan', yaxis='phase',
        spw='', field='', iteration='antenna', plotrange=[-1, -1, -180, 180])

plotcal(caltable='J0319+4130.D1', xaxis='antenna', yaxis='amp',
        figfile='plotcal_J0319+4130-D1.png')

# Solve for pol angle
polcal(vis=ms_active, caltable='3C138.X1',
       field='3C138', combine='scan',
       poltype='Xf', solint='inf',
       gaintable=["final_caltables/finaldelay.k",
                  "final_caltables/finalBPcal.b",
                  "J0319+4130_gaincal.g",
                  '3C138.Kcross',
                  'J0319+4130.D1'],
       gainfield=[])

plotcal(caltable='3C138.X1', xaxis='chan', yaxis='phase',
        figfile='plotcal_3c138-X1.png')

# Need to flux scale to J0319+4130. We have 2 standard cals here, so let's use 3C138,
# since we were already using it for the other polarization steps.
# myscale = fluxscale(vis=ms_active,
#                     caltable='3C138_gaincal.g',
#                     fluxtable='3C138.fluxscale',
#                     reference=['3C138'],
#                     transfer=['J0119+3210,J0319+4130'],
#                     incremental=False)

# Append on these three cal tables: Kcross, D1, X1
FinalGainTables = copy.copy(priorcals)
for i, field in enumerate(FinalGainTables):
    FinalGainTables[i] = os.path.join("final_caltables", field)
FinalGainTables.append('final_caltables/finaldelay.k')
FinalGainTables.append('final_caltables/finalBPcal.b')
# FinalGainTables.append('averagephasegain.g')
# FinalGainTables.append('finalampgaincal.g')
# FinalGainTables.append('finalphasegaincal.g')
# FinalGainTables.append('final_caltables/fluxgaincal.g')
# FinalGainTables.append("3C138.fluxscale")

default("applycal")

applycal(vis=ms_active,
         field='3C138',
         gaintable=FinalGainTables + ['3C138_gaincal.g',
                                      '3C138.Kcross',
                                      'J0319+4130.D1',
                                      '3C138.X1'],
         # [""] * (len(FinalGainTables) - 1) + ["3C138"] + [""] * 4,
         gainfield=[],
         # [""] * (len(FinalGainTables) - 1) + ["nearest"] + [""] * 4,
         interp=[],
         calwt=[False],
         parang=True,
         flagbackup=False)

applycal(vis=ms_active,
         field='J0319+4130',
         gaintable=FinalGainTables + ['J0319+4130_gaincal.g',
                                      '3C138.Kcross',
                                      'J0319+4130.D1',
                                      '3C138.X1'],
         # [""] * (len(FinalGainTables) - 1) + ["3C138"] + [""] * 4,
         gainfield=[],
         # [""] * (len(FinalGainTables) - 1) + ["nearest"] + [""] * 4,
         interp=[],
         calwt=[False],
         parang=True,
         flagbackup=False)


# Now use the gain tables from the pipeline for the non-pol cal sources.
FinalGainTables = copy.copy(priorcals)
for i, field in enumerate(FinalGainTables):
    FinalGainTables[i] = os.path.join("final_caltables", field)
FinalGainTables.append('final_caltables/finaldelay.k')
FinalGainTables.append('final_caltables/finalBPcal.b')
FinalGainTables.append('final_caltables/averagephasegain.g')
FinalGainTables.append('final_caltables/finalampgaincal.g')
FinalGainTables.append('final_caltables/finalphasegaincal.g')

applycal(vis=ms_active,
         field='J0119+3210',
         gaintable=FinalGainTables + ['3C138.Kcross',
                                      'J0319+4130.D1',
                                      '3C138.X1'],
         # [""] * (len(FinalGainTables) - 1) + ["J0119+3210"] + [""] * 4,
         gainfield=[],
         # [""] * (len(FinalGainTables) - 1) + ["nearest"] + [""] * 4,
         interp=[],
         calwt=[False],
         parang=True,
         flagbackup=False)

applycal(vis=ms_active,
         field='3C48',
         gaintable=FinalGainTables + ['3C138.Kcross',
                                      'J0319+4130.D1',
                                      '3C138.X1'],
         # [""] * (len(FinalGainTables) - 1) + ["J0119+3210"] + [""] * 4,
         gainfield=[],
         # [""] * (len(FinalGainTables) - 1) + ["nearest"] + [""] * 4,
         interp=[],
         calwt=[False],
         parang=True,
         flagbackup=False)

# Finally apply to the target
if "NGC604" in field_names:
    source_field = "NGC604"
elif "M33_Sarm" in field_names:
    source_field = "M33_Sarm"
else:
    raise ValueError("Cannot find a 'NGC604' or 'M33_Sarm' target field.")

applycal(vis=ms_active,
         field=source_field,
         gaintable=FinalGainTables +
         ['3C138.Kcross',
          'J0319+4130.D1',
          '3C138.X1'],
         # [""] * (len(FinalGainTables) - 1) + ["J0119+3210"] + [""] * 4,
         gainfield=[],
         # [""] * (len(FinalGainTables) - 1) + ["nearest"] + [""] * 4,
         interp=[],
         calwt=[False],
         parang=True,
         flagbackup=False)

# Now make some diagnostic plots
plotms(vis=ms_active, field='3C138', correlation='',
       antenna='', avgtime='60s',
       xaxis='channel', yaxis='amp', ydatacolumn='corrected',
       coloraxis='corr',
       plotfile='plotms_3C138-corrected-amp.png')

plotms(vis=ms_active, field='3C138', correlation='',
       antenna='', avgtime='60s',
       xaxis='channel', yaxis='phase', ydatacolumn='corrected',
       plotrange=[-1, -1, -180, 180], coloraxis='corr',
       plotfile='plotms_3C138-corrected-phase.png')

plotms(vis=ms_active, field='J0119+3210', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='channel', yaxis='amp', ydatacolumn='corrected',
       plotfile='plotms_J0119+3210-corrected-amp.png')

plotms(vis=ms_active, field='J0119+3210', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='channel', yaxis='phase', ydatacolumn='corrected',
       plotrange=[-1, -1, -180, 180], coloraxis='corr',
       plotfile='plotms_J0119+3210-corrected-phase.png')

plotms(vis=ms_active, field='J0319+4130', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='channel', yaxis='amp', ydatacolumn='corrected',
       plotfile='plotms_J0319+4130-corrected-amp.png')

plotms(vis=ms_active, field='J0319+4130', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='channel', yaxis='phase', ydatacolumn='corrected',
       plotrange=[-1, -1, -180, 180], coloraxis='corr',
       plotfile='plotms_J0319+4130-corrected-phase.png')

plotms(vis=ms_active, field='3C48', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='channel', yaxis='amp', ydatacolumn='corrected',
       plotfile='plotms_3C48-corrected-amp.png')

plotms(vis=ms_active, field='3C48', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='channel', yaxis='phase', ydatacolumn='corrected',
       plotrange=[-1, -1, -180, 180], coloraxis='corr',
       plotfile='plotms_3C48-corrected-phase.png')

# Amp vs phase

plotms(vis=ms_active, field='J0119+3210', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='phase', xdatacolumn='corrected', yaxis='amp',
       ydatacolumn='corrected',
       plotrange=[], coloraxis='corr',
       plotfile='plotms_J0119+3210-corrected-ampvsphase.png')

plotms(vis=ms_active, field='J0319+4130', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='phase', xdatacolumn='corrected', yaxis='amp',
       ydatacolumn='corrected',
       plotrange=[], coloraxis='corr',
       plotfile='plotms_J0319+4130-corrected-ampvsphase.png')

plotms(vis=ms_active, field='3C48', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='phase', xdatacolumn='corrected', yaxis='amp',
       ydatacolumn='corrected',
       plotrange=[], coloraxis='corr',
       plotfile='plotms_3C48-corrected-ampvsphase.png')

plotms(vis=ms_active, field='3C138', correlation='RR,LL',
       timerange='', antenna='', avgtime='60s',
       xaxis='phase', xdatacolumn='corrected', yaxis='amp',
       ydatacolumn='corrected',
       plotrange=[], coloraxis='corr',
       plotfile='plotms_3C138-corrected-ampvsphase.png')
