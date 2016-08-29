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

'''
Facilitates the restart of the pipeline.
EVLA_pipe_restore.py needs to be run before hand

'''

pipeline_scripts = ['startup', 'import', 'hanning', 'msinfo', 'flagall',
                    'calprep', 'priorcals', 'testBPdcals',
                    'flag_baddeformatters', 'flag_baddeformattersphase',
                    'checkflag', 'semiFinalBPdcals1',
                    'checkflag_semiFinal', 'semiFinalBPdcals2',
                    'solint', 'testgains', 'fluxgains', 'fluxboot',
                    'finalcals', 'applycals', 'targetflag',
                    'statwt', 'plotsummary', 'filecollect', 'weblog']

#last script that was started
last_state = time_list[-1]['pipestate']
last_status = time_list[-1]['status']

#Find index where we need to pick up at
script_index = [i for i, x in enumerate(pipeline_scripts) if last_state == x][0]

#If script ended successfully then start on next script
if (last_status == 'end'):
    script_index = script_index + 1

if (last_status == 'start'):
    time_list.pop(-1)

for i in range( script_index, len(pipeline_scripts) ):
    try:
        execfile(pipepath+'EVLA_pipe_'+pipeline_scripts[i]+'.py')
    except Exception, e:
        logprint ("Exiting script: "+str(e))
    except KeyboardInterrupt, e:
        logprint ("Keyboard Interrupt: "+str(e))
