######################################################################
#
# Copyright (C) 2013
# Associated Universities, Inc. Washington DC, USA,
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning VLA Pipelines should be addressed as follows:
#    Please register and submit helpdesk tickets via: https://help.nrao.edu
#    Postal address:
#              National Radio Astronomy Observatory
#              VLA Pipeline Support Office
#              PO Box O
#              Socorro, NM,  USA
#
######################################################################

# PRIOR CALIBRATIONS

# Gain curves, opacities, antenna position corrections
# requantizer gains (will require CASA 4.1)

logprint ("Starting EVLA_pipe_priorcals.py", logfileout='logs/priorcals.log')
time_list=runtiming('priorcals', 'start')
QA2_priorcals='Pass'

# Gain curves first:

default(gencal)
vis=ms_active
caltable='gain_curves.g'
caltype='gc'
spw=''
antenna=''
pol=''
parameter=[]
gencal()

# Opacities:

default(gencal)
vis=ms_active
caltable='opacities.g'
caltype='opac'
spw=all_spw
antenna=''
pol=''
parameter=tau
gencal()

# Apply switched power calibration (when commissioned); for now, just
# requantizer gains (needs casa4.1!), and only for data with
# sensible switched power tables (Feb 24, 2011)

if startdate >= 55616.6:
    default(gencal)
    vis=ms_active
    caltable='requantizergains.g'
    caltype='rq'
    spw=''
    antenna=''
    pol=''
    parameter=[]
    gencal()

if os.path.exists('requantizergains.g'):
    priorcals=['gain_curves.g','opacities.g','requantizergains.g']
else:
    priorcals=['gain_curves.g','opacities.g']

# Correct for antenna position errors, if known

# NB: for the realtime pipeline these will not be available yet, but all
# SBs that do not have good antenna positions should be re-processed when
# they are available

try:
    default(gencal)
    vis=ms_active
    caltable='antposcal.p'
    caltype='antpos'
    spw=''
    antenna=''
    pol=''
    parameter=[]
    gencal()
#
    if os.path.exists('antposcal.p'):
        priorcals.append('antposcal.p')
        antenna_offsets=correct_ant_posns(ms_active)
        logprint ("Correcting for known antenna position errors", logfileout='logs/priorcals.log')
        logprint (str(antenna_offsets), logfileout='logs/priorcals.log')
    else:
        logprint ("No antenna position corrections found/needed", logfileout='logs/priorcals.log')
except:
    logprint ("No antenna position corrections found/needed", logfileout='logs/priorcals.log')

# Lastly, make switched power table.  This is not used in the
# pipeline, but may be used for QA and for flagging, especially at
# S-band for fields near the geostationary satellite belt.  Only
# relevant for data taken on 24-Feb-2011 or later.

if startdate >= 55616.6:
    default(gencal)
    vis=ms_active
    caltable='switched_power.g'
    caltype='swpow'
    spw=''
    antenna=''
    pol=''
    parameter=[]
    gencal()
#
# plot switched power table
#
    logprint ("Plotting switched power table", logfileout='logs/priorcals.log')
#
    nplots=int(numAntenna/3)
#
    if ((numAntenna%3)>0):
        nplots = nplots + 1

    for ii in range(nplots):
        filename='switched_power'+str(ii)+'.png'
        syscommand='rm -rf '+filename
        os.system(syscommand)
#
        antPlot=str(ii*3)+'~'+str(ii*3+2)
#
        default('plotcal')
        caltable='switched_power.g'
        xaxis='time'
        yaxis='spgain'
        poln=''
        field=''
        antenna=antPlot
        spw=''
        timerange=''
        subplot=311
        overplot=False
        clearpanel='Auto'
        iteration='antenna'
        plotrange=[]
        showflags=False
        plotsymbol='o'
        plotcolor='blue'
        markersize=5.0
        fontsize=10.0
        showgui=False
        figfile=filename
        plotcal()

    for ii in range(nplots):
        filename='Tsys'+str(ii)+'.png'
        syscommand='rm -rf '+filename
        os.system(syscommand)
#
        antPlot=str(ii*3)+'~'+str(ii*3+2)
#
        default('plotcal')
        caltable='switched_power.g'
        xaxis='time'
        yaxis='tsys'
        poln=''
        field=''
        antenna=antPlot
        spw=''
        timerange=''
        subplot=311
        overplot=False
        clearpanel='Auto'
        iteration='antenna'
        plotrange=[]
        showflags=False
        plotsymbol='o'
        plotcolor='blue'
        markersize=5.0
        fontsize=10.0
        showgui=False
        figfile=filename
        plotcal()

# Until we know what error messages to search for in priorcals,
# leave QA2 score set to "Pass".

logprint ("QA2 score: "+QA2_priorcals, logfileout='logs/priorcals.log')
logprint ("Finished EVLA_pipe_priorcals.py", logfileout='logs/priorcals.log')
time_list=runtiming('priorcals', 'end')

pipeline_save()
