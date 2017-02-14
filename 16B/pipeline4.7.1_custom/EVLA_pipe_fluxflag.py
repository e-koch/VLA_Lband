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

# FLAG GAIN TABLE PRIOR TO FLUX DENSITY BOOTSTRAPPING

# Flag fluxgaincal.g interactively to remove bad times/antennas
# NB: not part of the regular pipeline script, must not do runtiming
# on it otherwise EVLA_pipe_restart.py doesn't work

logprint ("Starting EVLA_pipe_fluxflag.py", logfileout='logs/fluxflag.log')
#time_list=runtiming('fluxflag', 'start')

default('plotcal')
caltable='fluxgaincal.g'
xaxis='time'
yaxis='amp'
subplot=111
iteration=''
plotcal()


logprint ("Finished EVLA_pipe_fluxflag.py", logfileout='logs/fluxflag.log')
#time_list=runtiming('fluxflag', 'end')

pipeline_save()
