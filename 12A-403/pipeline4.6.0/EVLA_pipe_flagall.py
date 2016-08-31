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

# FLAG SOME STUFF WE KNOW ABOUT (TIME-BASED)
# online flags, zeroes, pointing scans, quacking
# Use flagdata in list mode
# Do zero and shadow flagging separately STM 20121108
#
# Note: Urvashi Rao points out that it may be possible to
# consolidate these to run faster (minus the flagcmd step) using the
# toolkit, see
# http://www.aoc.nrao.edu/~rurvashi/ActiveFlaggerDocs/node11.html

# Apply deterministic flags

logprint ("Starting EVLA_pipe_flagall.py", logfileout='logs/flagall.log')
time_list=runtiming('flagall', 'start')
QA2_flagall='Pass'

logprint ("Deterministic flagging", logfileout='logs/flagall.log')

outputflagfile = 'flagging_commands1.txt'
syscommand='rm -rf '+outputflagfile
os.system(syscommand)

logprint ("Determine fraction of time on-source (may not be correct for pipeline re-runs on datasets already flagged)", logfileout='logs/flagall.log')

# report initial statistics
default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
myinitialflags = flagdata()
#clearstat()
logprint ("Initial flags summary", logfileout='logs/flagall.log')

start_total = myinitialflags['total']
start_flagged = myinitialflags['flagged']
logprint ("Initial flagged fraction = "+str(start_flagged/start_total), logfileout='logs/flagall.log')

default('flagcmd')
vis=ms_active
inpmode='xml'
tbuff=1.5*int_time
ants=''
reason='ANTENNA_NOT_ON_SOURCE'
action='apply'
flagbackup=True
savepars=True
outfile=outputflagfile
flagcmd()
logprint ("ANTENNA_NOT_ON_SOURCE flags carried out", logfileout='logs/flagall.log')

# Now shadow flagging
default('flagdata')
vis=ms_active
mode='shadow'
tolerance=0.0
action='apply'
flagbackup=False
savepars=False
flagdata()
#clearstat()

# report new statistics
default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
slewshadowflags=flagdata()

init_on_source_vis = start_total-slewshadowflags['flagged']

logprint ("Initial on-source fraction = "+str(init_on_source_vis/start_total), logfileout='logs/flagall.log')

# Restore original flags

default('flagmanager')
vis=ms_active
mode='restore'
versionname='flagcmd_1'
merge='replace'
flagmanager()

syscommand='rm -rf '+outputflagfile
os.system(syscommand)

# First do zero flagging (reason='CLIP_ZERO_ALL')
default('flagdata')
vis=ms_active
mode='clip'
correlation='ABS_ALL'
clipzeros=True
action='apply'
flagbackup=False
savepars=False
outfile=outputflagfile
myzeroflags = flagdata()
#clearstat()
logprint ("Zero flags carried out", logfileout='logs/flagall.log')

# Now report statistics
default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
myafterzeroflags = flagdata()
#clearstat()
logprint ("Zero flags summary", logfileout='logs/flagall.log')

afterzero_total = myafterzeroflags['total']
afterzero_flagged = myafterzeroflags['flagged']
logprint ("After ZERO flagged fraction = "+str(afterzero_flagged/afterzero_total), logfileout='logs/flagall.log')

zero_flagged = myafterzeroflags['flagged'] - myinitialflags['flagged']
logprint ("Delta ZERO flagged fraction = "+str(zero_flagged/afterzero_total), logfileout='logs/flagall.log')

# Now shadow flagging
default('flagdata')
vis=ms_active
mode='shadow'
tolerance=0.0
action='apply'
flagbackup=False
savepars=False
flagdata()
#clearstat()
logprint ("Shadow flags carried out", logfileout='logs/flagall.log')

# Now report statistics after shadow
default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
myaftershadowflags = flagdata()
#clearstat()
logprint ("Shadow flags summary", logfileout='logs/flagall.log')

aftershadow_total = myaftershadowflags['total']
aftershadow_flagged = myaftershadowflags['flagged']
logprint ("After SHADOW flagged fraction = "+str(aftershadow_flagged/aftershadow_total), logfileout='logs/flagall.log')

shadow_flagged = myaftershadowflags['flagged'] - myafterzeroflags['flagged']
logprint ("Delta SHADOW flagged fraction = "+str(shadow_flagged/aftershadow_total), logfileout='logs/flagall.log')

default('flagcmd')
vis=ms_active
inpmode='xml'
tbuff=1.5*int_time
ants=''
reason='any'
action='apply'
flagbackup=False
savepars=True
outfile=outputflagfile
flagcmd()
logprint ("Online flags applied", logfileout='logs/flagall.log')

#Define list of flagdata parameters to use in 'list' mode
flagdata_list=[]
cmdreason_list=[]

# Flag pointing scans, if there are any
if (len(pointing_state_IDs) != 0):
    logprint ("Flag pointing scans", logfileout='logs/flagall.log')
    flagdata_list.append("mode='manual' intent='*POINTING*'")
    cmdreason_list.append('pointing')
   
# Flag setup scans
logprint ("Flag setup scans", logfileout='logs/flagall.log')
flagdata_list.append("mode='manual' intent='UNSPECIFIED#UNSPECIFIED'")
cmdreason_list.append('setup')

logprint ("Flag setup scans", logfileout='logs/flagall.log')
flagdata_list.append("mode='manual' intent='SYSTEM_CONFIGURATION#UNSPECIFIED'")
cmdreason_list.append('setup')

# Quack the data
logprint ("Quack the data", logfileout='logs/flagall.log')
flagdata_list.append("mode='quack' scan=" + quack_scan_string +
    " quackinterval=" + str(1.5*int_time) + " quackmode='beg' " +
    "quackincrement=False")
cmdreason_list.append('quack')

######################################################################

# FLAG SOME MORE STUFF (CHANNEL-BASED)
# Flag end 3 channels of each spw
logprint ("Flag end 5 percent of each spw or minimum of 3 channels", logfileout='logs/flagall.log')

SPWtoflag=''

for ispw in range(numSpws):
    fivepctch=int(0.05*channels[ispw])
    startch1=0
    startch2=fivepctch-1
    endch1=channels[ispw]-fivepctch
    endch2=channels[ispw]-1
    
    #Minimum number of channels flagged must be three on each end
    if (fivepctch < 3):
        startch2=2
        endch1=channels[ispw]-3
    
    if (ispw<max(range(numSpws))):
        SPWtoflag=SPWtoflag+str(ispw)+':'+str(startch1)+'~'+str(startch2)+';'+str(endch1)+'~'+str(endch2)+','
    else:
        SPWtoflag=SPWtoflag+str(ispw)+':'+str(startch1)+'~'+str(startch2)+';'+str(endch1)+'~'+str(endch2)

flagdata_list.append("mode='manual' spw='"+SPWtoflag+"'")
cmdreason_list.append('spw_ends')

# Flag 10 end channels at edges of basebands
#
# NB: assumes continuum set-up that fills baseband; will want to modify
# for narrow spws or spectroscopy!
#
bottomSPW=''
topSPW=''

for ii in range(0,len(low_spws)):
    if (ii == 0):
        bspw=low_spws[ii]
        tspw=high_spws[ii]
        endch1=channels[tspw]-10
        endch2=channels[tspw]-1
        bottomSPW=str(bspw)+':0~9'
        topSPW=str(tspw)+':'+str(endch1)+'~'+str(endch2)
    else:
        bspw=low_spws[ii]
        tspw=high_spws[ii]
        endch1=channels[tspw]-10
        endch2=channels[tspw]-1
        bottomSPW=bottomSPW+','+str(bspw)+':0~9'
        topSPW=topSPW+','+str(tspw)+':'+str(endch1)+'~'+str(endch2)
 
if (bottomSPW != ''):
    logprint ("Flag end 10 channels at edges of basebands", logfileout='logs/flagall.log')
    SPWtoflag=bottomSPW+','+topSPW
    flagdata_list.append("mode='manual' spw='"+SPWtoflag+"'")
    cmdreason_list.append('baseband_edge_chans')

######################################################################

#Write out list for use in flagdata mode 'list'
f = open(outputflagfile, 'a')
for line in flagdata_list:
    f.write(line+"\n")
f.close()

# Apply all flags
logprint ("Applying all flags to data", logfileout='logs/flagall.log')

default('flagdata')
vis=ms_active
mode='list'
inpfile=outputflagfile
action='apply'
flagbackup=False
savepars=True
cmdreason=string.join(cmdreason_list, ',')
flagdata()
#clearstat()

logprint ("Flagging completed ", logfileout='logs/flagall.log')
logprint ("Flag commands saved in file "+outputflagfile, logfileout='logs/flagall.log')

# Save flags
logprint ("Saving flags", logfileout='logs/flagall.log')


default('flagmanager')
vis=ms_active
mode='save'
versionname='allflags1'
comment='Deterministic flags saved after application'
merge='replace'
flagmanager()
logprint ("Flag column saved to "+versionname, logfileout='logs/flagall.log')


# report new statistics
default('flagdata')
vis=ms_active
mode='summary'
spwchan=True
spwcorr=True
basecnt=True
action='calculate'
savepars=False
all_flags = flagdata()

frac_flagged_on_source1 = 1.0-((start_total-all_flags['flagged'])/init_on_source_vis)

logprint ("Fraction of on-source data flagged = "+str(frac_flagged_on_source1), logfileout='logs/flagall.log')

if (frac_flagged_on_source1 >= 0.3):
    QA2_flagall='Fail'

logprint ("Finished EVLA_pipe_flagall.py", logfileout='logs/flagall.log')
logprint ("QA2 score: "+QA2_flagall, logfileout='logs/flagall.log')
time_list=runtiming('flagall', 'end')

pipeline_save()
