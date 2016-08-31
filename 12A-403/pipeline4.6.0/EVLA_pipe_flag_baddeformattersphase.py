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

# Uncomment the following to run without rest of pipeline
#def logprint(printstr,logfileout):
#    print printstr
#    return
#
#from lib_EVLApipeutils import getBCalStatistics

######################################################################

# Determine bad deformatters in the MS and flag them
# Looks for bandpass solutions that have small ratio of min/max amplitudes
# Version updated 2012-12-05 STM
# Version updated 2012-12-07 STM - check number of spw per baseband
# Version updated 2012-12-17 STM - add phase checks: requires 20121217 lib_EVLApipeutils.py

#!cp /home/sandrock2/smyers/casa/pipeline/pipeline4.1/EVLA_pipe_flag_baddeformattersphase.py .

logprint("Starting EVLA_pipe_flag_baddeformattersphase.py", logfileout='logs/flag_baddeformattersphase.log')
time_list=runtiming('flag_baddeformattersphase', 'start')
QA2_flag_baddeformattersphase='Pass'

# Set control parameters here
# Print detailed flagging stats
doprintall = True
# Which quantity to test? ['amp','phase','real','imag']
testq = 'phase'
# Which stat to use? ['min','max','mean','var'] or 'rat'=min/max or 'diff'=max-min
tstat = 'diff'
# Limit for test (flag values under/over this limit)
testlimit = 50
testunder = False
# Number of spw per baseband to trigger flagging entire baseband
nspwlimit = 4
# Flag individual spws when below nspwlimit
doflagundernspwlimit = True
# Flag data for spws with no unflagged channel solutions in any poln?
doflagemptyspws = False
# Actually flag the data based on the derived flags (or just report)?

if startdate <= 56062.7:
    doflagdata = False
else:
    doflagdata = True

# Define the table to run this on
calBPtablename ='testBPcal.b'
# Define the REASON given for the flags
flagreason = 'bad_deformatters_phase or RFI'

logprint("Will test on quantity: "+testq, logfileout='logs/flag_baddeformattersphase.log')
logprint("Will test using statistic: "+tstat, logfileout='logs/flag_baddeformattersphase.log')
if testunder:
    logprint("Will flag values under limit = "+str(testlimit), logfileout='logs/flag_baddeformattersphase.log')
else:
    logprint("Will flag values over limit = "+str(testlimit), logfileout='logs/flag_baddeformattersphase.log')

logprint("Will identify basebands with more than "+str(nspwlimit)+" bad spw", logfileout='logs/flag_baddeformattersphase.log')

if doflagundernspwlimit:
    logprint("Will identify individual spw when less than "+str(nspwlimit)+" bad spw", logfileout='logs/flag_baddeformattersphase.log')

if doflagemptyspws:
    logprint("Will identify spw with no unflagged channels", logfileout='logs/flag_baddeformattersphase.log')

logprint("Will use flag REASON = "+flagreason, logfileout='logs/flag_baddeformattersphase.log')

if doflagdata:
    logprint("Will flag data based on what we found", logfileout='logs/flag_baddeformattersphase.log')
else:
    logprint("Will NOT flag data based on what we found", logfileout='logs/flag_baddeformattersphase.log')

calBPstatresult = getBCalStatistics(calBPtablename)
flaglist = []
extflaglist = []
for iant in calBPstatresult['antband'].keys():
    antName = calBPstatresult['antDict'][iant]
    badspwlist = []
    flaggedspwlist = []
    for rrx in calBPstatresult['antband'][iant].keys():
        for bband in calBPstatresult['antband'][iant][rrx].keys():
            #
            # List of spw in this baseband
            spwl = calBPstatresult['rxBasebandDict'][rrx][bband]
            #
            nbadspws = 0
            badspws = []
            flaggedspws = []
            if len(spwl)>0:
                if doprintall:
                    strprt=' Ant %s (%s) %s %s processing spws=%s' % (str(iant),antName,rrx,bband,str(spwl))
                    logprint (strprt, logfileout='logs/flag_baddeformattersphase.log')
                for ispw in spwl:
                    #
                    testvalid = False
                    if calBPstatresult['antspw'][iant].has_key(ispw):
                        for poln in calBPstatresult['antspw'][iant][ispw].keys():
                            # Get stats of this ant/spw/poln
                            nbp = calBPstatresult['antspw'][iant][ispw][poln]['inner']['number']
                            #
                            if nbp>0:
                                if tstat=='rat':
                                    bpmax = calBPstatresult['antspw'][iant][ispw][poln]['inner'][testq]['max']
                                    bpmin = calBPstatresult['antspw'][iant][ispw][poln]['inner'][testq]['min']
                                    #
                                    if bpmax==0.0:
                                        tval = 0.0
                                    else:
                                        tval = bpmin/bpmax
                                elif tstat=='diff':
                                    bpmax = calBPstatresult['antspw'][iant][ispw][poln]['inner'][testq]['max']
                                    bpmin = calBPstatresult['antspw'][iant][ispw][poln]['inner'][testq]['min']
                                    #
                                    tval = bpmax-bpmin
                                else:
                                    # simple test on quantity
                                    tval = calBPstatresult['antspw'][iant][ispw][poln]['inner'][testq][tstat]
                                if not testvalid:
                                    testval = tval
                                    testvalid = True
                                elif testunder:
                                    if tval<testval:
                                        testval = tval
                                else:
                                    if tval>testval:
                                        testval = tval
                                #
                        # Test on extrema of the polarizations for this ant/spw
                        if not testvalid:
                            # these have no unflagged channels in any poln
                            flaggedspws.append(ispw)
                        else:
                            if (testunder and testval<testlimit) or (not testunder and testval>testlimit):
                                nbadspws += 1
                                badspws.append(ispw)
                                if doprintall:
                                    strprt='  Found Ant %s (%s) %s %s spw=%s %s %s=%6.4f' % (str(iant),antName,rrx,bband,str(ispw),testq,tstat,testval)
                                    logprint (strprt, logfileout='logs/flag_baddeformattersphase.log')
                    else:
                        # this spw is missing from this antenna/rx
                        if doprintall:
                            strprt='  Ant %s (%s) %s %s spw=%s missing solution' % (str(iant),antName,rrx,bband,str(ispw))
                            logprint (strprt, logfileout='logs/flag_baddeformattersphase.log')
                #
            #
            # Test to see if this baseband should be entirely flagged
            if nbadspws>0 and nbadspws>=nspwlimit:
                # Flag all spw in this baseband
                bbspws = calBPstatresult['rxBasebandDict'][rrx][bband]
                badspwlist.extend(bbspws)
                strprt='Ant %s (%s) %s %s bad baseband spws=%s' % (str(iant),antName,rrx,bband,str(bbspws))
                logprint (strprt, logfileout='logs/flag_baddeformattersphase.log')
            elif nbadspws>0 and doflagundernspwlimit:
                # Flag spws individually
                badspwlist.extend(badspws)
                strprt='Ant %s (%s) %s %s bad spws=%s' % (str(iant),antName,rrx,bband,str(badspws))
                logprint (strprt, logfileout='logs/flag_baddeformattersphase.log')
                print ""
            if len(flaggedspws)>0:
                # these spws have no unflagged channels in any pol
                flaggedspwlist.extend(flaggedspws)
                strprt='Ant %s (%s) %s %s no unflagged solutions spws=%s ' % (str(iant),antName,rrx,bband,str(flaggedspws))
                logprint (strprt, logfileout='logs/flag_baddeformattersphase.log')
                #
            #
        #
    #
    if len(badspwlist)>0:
        spwstr = '' 
        for ispw in badspwlist:
            if spwstr=='':
                spwstr = str(ispw)
            else:
                spwstr+=','+str(ispw)
        #
        #reastr = 'bad_deformatters'
        reastr = flagreason
        # Add entry for this antenna
        #flagstr = "antenna='"+str(iant)+"' spw='"+spwstr+"' reason='"+reastr+"'"
        # Use name for flagging
        flagstr = "mode='manual' antenna='"+antName+"' spw='"+spwstr+"'"
        flaglist.append(flagstr)

    if doflagemptyspws and len(flaggedspwlist)>0:
        spwstr = '' 
        for ispw in flaggedspwlist:
            if spwstr=='':
                spwstr = str(ispw)
            else:
                spwstr+=','+str(ispw)
        #
        # Add entry for this antenna
        reastr = 'no_unflagged_solutions'
        #flagstr = "antenna='"+str(iant)+"' spw='"+spwstr+"' reason='"+reastr+"'"
        # Use name for flagging
        flagstr = "mode='manual' antenna='"+antName+"' spw='"+spwstr+"'"
        extflaglist.append(flagstr)

nflagcmds = len(flaglist)+len(extflaglist)
if nflagcmds<1:
    logprint ("No bad basebands/spws found", logfileout='logs/flag_baddeformattersphase.log')
else:
    logprint ("Possible bad basebands/spws found:", logfileout='logs/flag_baddeformattersphase.log')
    #
    for flagstr in flaglist:
        logprint ("   "+flagstr, logfileout='logs/flag_baddeformattersphase.log')
    if len(extflaglist)>0:
        logprint ("   ", logfileout='logs/flag_baddeformattersphase.log')
        for flagstr in extflaglist:
            logprint ("   "+flagstr, logfileout='logs/flag_baddeformattersphase.log')
        flaglist.extend(extflaglist)
    if doflagdata:
        logprint ("Flagging these in the ms:", logfileout='logs/flag_baddeformattersphase.log')
        #Now flag with flaglist
        default('flagdata')
        vis=ms_active
        mode='list'
        inpfile=flaglist
        action='apply'
        flagbackup=True
        savepars=True
        flagdata()

# Until we know what the QA criteria are for this script, leave QA2
# set score to "Pass".
 
logprint ("Finished EVLA_pipe_flag_baddeformattersphase.py", logfileout='logs/flag_baddeformattersphase.log')
time_list=runtiming('flag_baddeformattersphase', 'end')

pipeline_save()
