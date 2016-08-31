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
Script utility for moving relevant plots and tables after calibration 
    is complete
VLA pipeline
August 2012, B. Kent, NRAO
November 2012, B. Butler, NRAO
'''
import os
import stat
import glob
import shutil
import copy

######################################################################

logprint ('Starting EVLA_pipe_filecollect.py', logfileout='logs/filecollect.log')
time_list=runtiming('filecollect', 'start')

caltables_dir = './final_caltables'
weblog_dir = './weblog'

#Check if directories exists - if not create one
if not os.path.exists(caltables_dir):
    os.makedirs(caltables_dir)

if not os.path.exists(weblog_dir):
    os.makedirs(weblog_dir)

#Move all png files/plots
png_files = glob.glob('./*.png')
for file in png_files:
    try:
        shutil.move(file, weblog_dir+'/.')
    except:
        logprint('Unable to move ' + file, logfileout='logs/filecollect.log')

logprint('Plots moved to ' + weblog_dir, logfileout='logs/filecollect.log')

#Listobs output
listobs_output = glob.glob('./*.listobs')
for file in listobs_output:
    try:
        shutil.move(file, weblog_dir+'/.')
    except:
        logprint('Unable to move ' + file, logfileout='logs/filecollect.log')
  
#Move calibration tables into caltables_dir
cal_files=copy.copy(priorcals)
cal_files.append('switched_power.g')
cal_files.append('fluxgaincal.g')
cal_files.append('finaldelay.k')
cal_files.append('finalBPcal.b')
cal_files.append('averagephasegain.g')
cal_files.append('finalampgaincal.g')
cal_files.append('finalphasegaincal.g')

for file in cal_files:
    try:
        shutil.move(file, caltables_dir+'/.')
    except:
        logprint('Unable to move '+file,logfileout='logs/filecollect.log')
  
logprint('Final calibration tables moved to ' + caltables_dir, logfileout='logs/filecollect.log')

######################################################################
#Pickle up the timing profile list
try:
    gmt_time = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    file_time = strftime("%d%b%Y_%H%M%Sgmt", gmtime())
    
    #compute size of ms directory
    ms_size = 0
    try:
        for (path, dirs, files) in os.walk(ms_active):
            for file in files:
                filename = os.path.join(path, file)
                ms_size += os.path.getsize(filename)
        #Size of the ms in Gigabytes
        ms_size = ms_size / (1024.0*1024.0*1024.0) 
    except:
        logprint('Unable to determine size of ms on disk')
    
    pipeprofile = {'SDM_name':SDM_name, 
               'time_list':time_list, 
               'gmt_time':gmt_time,
               'version':version,
               'svnrevision':svnrevision,
               'ms_size':ms_size}
               
    if os.path.exists('logs'):
        logprint('logs directory exists')
    else:
        logprint('logs directory does not exist / not in path')
    
#
# if SDM_name is a full path (as it is for automatic pipeline
# executions, it will have '/'s in it.  search for the right-most one,
# and only use the part after that if it exists.  otherwise, just use
# the full name (since it's one done by hand).
#
    right_index = SDM_name.rfind('/')
    pickle_filename = 'logs/profile_' + SDM_name[right_index+1:] + '_' + file_time + '.p'
    logprint ('writing pickle file: ' + pickle_filename, logfileout='logs/completion.log')
    try:         
        istat = pickle.dump(pipeprofile, open(pickle_filename, 'wb'))
        logprint('Pickle dump of profile successful')
    except:
        logprint('Unable to dump pickle file')
    #Load this dict with pipeprofile = pickle.load( open( "<file.p>", "rb" ) )

    logprint ('Timing profile written to logs/timing.log', logfileout='logs/completion.log')
    logprint ('Completed on ' + gmt_time + ' with pipeline SVN revision ' + svnrevision, logfileout='logs/timing.log')
    
except:
    logprint ('Problem writing timing profile')
    
    
######################################################################
#Attempt to copy files to /lustre/aoc/cluster/pipeline/stats
stats_dir = '/lustre/aoc/cluster/pipeline/stats'

if os.path.exists(stats_dir):
    log_filename = stats_dir + '/profile_' + SDM_name[right_index+1:] + '_' + file_time + '.log'
    try:
        shutil.copyfile('logs/timing.log', log_filename)
        logprint ('Copied timing log to ' + stats_dir)
    except:
        logprint('Unable to copy timing log to ' + stats_dir)

    profile_filename = stats_dir + '/profile_' + SDM_name[right_index+1:] + '_' + file_time+ '.p'
    try:
        shutil.copyfile(pickle_filename, profile_filename)
        logprint ('Copied profile to ' + stats_dir)
    except:
        logprint('Unable to copy profile to ' + stats_dir)
else:
    logprint (stats_dir + ' does not exist or not accessible')

#
# if this is an automatic execution, try to fix file permissions, so 
# folks in the vlapipe group can clean up after things are moved into
# the Trash folder...
#

if SDM_name_already_defined:
#
# if it turns out things in the rawdata or results directories are
# improperly permissioned, uncomment the getcwd() and chdir() foo
#
#   cwd = os.getcwd()
#   os.chdir('..')
    for path, dirs, files in os.walk('.'):
        for dir in dirs:
            full_dir_name = os.path.join(path, dir)
            st = os.stat(full_dir_name)
            if not (bool(st.st_mode & stat.S_IRGRP) and bool(st.st_mode & stat.S_IWGRP) and bool(st.st_mode & stat.S_IXGRP)):
                os.system("chmod -f g+rwx %s" % full_dir_name)
        for file in files:
            if file != 'epilogue.sh':
                full_file_name = os.path.join(path, file)
                st = os.stat(full_file_name)
                if not (bool(st.st_mode & stat.S_IRGRP) and bool(st.st_mode & stat.S_IWGRP) and bool(st.st_mode & stat.S_IXGRP)):
                    os.system("chmod -f g+rwx %s" % full_file_name)
#   os.chdir(cwd)

logprint ('Finished EVLA_pipe_filecollect.py', logfileout='logs/filecollect.log')
time_list=runtiming('filecollect', 'end')

pipeline_save()




