######################################################################
# Functions for EVLA pipeline script
# Runs on CASA 4.1.0
# 02/05/13 C. Chandler
# 09/23/14 C. Chandler (updated for CASA 4.2.2)
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

# From ALMA/Nick Elias:

# Usage:
# out = selectReferenceAntenna( '<ms>', '<spws>', '<fields>' )

# ------------------------------------------------------------------------------

# selectReferenceAntenna

# Description:
# ------------
# This function selects the reference antenna based on x-y position within the
# the array as well as the number of completely flagged rows.

# NB: This function is meant to be a temporary stop gap.  The pipeline version
# will be organized differently and perhaps include additional criteria.

# NB: This function determines the reference antenna by separate distance and
# flag "scores" and creating from them a combined "score."

# NB: This function does not use intents to determine which fields are
# calibrators.  The user must select fields.  

# NB: This function does look for phase jumps.

# Inputs:
# -------
# ms     - The MS name, specified as a python string.  If the MS is not in the
#          present working directory, the relative or absolute paths must be
#          specified.
# spws   - The spectral window IDs, specified as a comma-delimited python string.
#          The default is '' (all spectral window IDs).  NB: Do not use names
#          from the NAME column in the SPECTRAL_WINDOW subtable, use only the
#          numerical IDs.  Also, do no include channel selection.
# fields - The field IDs, specified as a comma-delimited python string.  The
#          default is '' (all field IDs).  NB: Do not use names from the NAME
#          column in the FIELD subtable.  Use only the numerical IDs.

# Outputs:
# --------
# A python combination, returned via the function value, whose elements are:
# out[0] - A python string containing the name of the reference antenna ID (not
#          the numerical antenna IDs).
# out[1] - A numpy array of strings containing the antenna IDs, ordered
#          according to total score (the first element is the reference antenna
#          ID).
# out[2] - A python dictionary containing the total scores.  The keywords are the
#          antenna IDs and the elements are the total scores.
# out[3] - A python dictionary containing the distance scores.  The keywords are
#          the antenna IDs and the elements are the distance scores.
# out[4] - A python dictionary containing the flag row scores.  The keywords are
#          the antenna IDs and the elements are the flag row scores.

# Modification history:
# ---------------------
# 2012 Feb 13 - Nick Elias, NRAO
#               This function was pilfered from the original pipeline code and
#               cleaned up.  It was converted from a member function to a global
#               function.  All references to member functions were removed.  The
#               spectral window and field inputs are now comma-delimited strings
#               with '' defaults (before, they were lists of integers with no
#               defaults).

# ------------------------------------------------------------------------------

def selectReferenceAntenna( ms, spws='', fields='' ):

	# Create the local instance of the table, measures, and quanta tools

	tbLoc = casa.__tablehome__.create()
	meLoc = casa.__measureshome__.create()
	qaLoc = casa.__quantahome__.create()


	# Fix the spectral window input

	if spws != '':

		spwSList = spws.split( ',' )

		spwList = []
		for s in spwSList: spwList.append( int(s) )

	else:

		tbLoc.open( ms + '/SPECTRAL_WINDOW' )
		rSPW = range( tbLoc.nrows() )
		tbLoc.close()

		spwList = []
		for s in rSPW: spwList.append( int(s) )


	# Fix the field input

	if fields != '':

		fieldSList = fields.split( ',' )

		fieldList = []
		for f in fieldSList: fieldList.append( int(f) )

	else:

		tbLoc.open( ms + '/FIELD' )
		rField = range( tbLoc.nrows() )
		tbLoc.close()

		fieldList = []
		for f in rField: fieldList.append( int(f) )


	# Read the antenna positions.  These seem to be in metres rel to centre
	# of the Earth.
     
	tbLoc.open( '%s/ANTENNA' % ms )

	position = tbLoc.getcol( 'POSITION' )
	flag_row = tbLoc.getcol( 'FLAG_ROW' )
	name = tbLoc.getcol( 'NAME' )
	position_keywords = tbLoc.getcolkeywords( 'POSITION' )

	tbLoc.close()


	# Make a position Measure for each antenna, this stores info in
	# terms of long, lat and distance from centre of Earth.

	antenna_position = {}

	for row,ant in enumerate( name ):
		if not flag_row[row]:
			rf = position_keywords['MEASINFO']['Ref']
			v0 = qaLoc.quantity( position[0,row],
			    position_keywords['QuantumUnits'][0] )
			v1 = qaLoc.quantity( position[1,row],
			    position_keywords['QuantumUnits'][1] )
			v2 = qaLoc.quantity( position[2,row],
			    position_keywords['QuantumUnits'][2] )
			ant_pos = meLoc.position( rf=rf, v0=v0, v1=v1, v2=v2 )
			antenna_position[ant] = ant_pos


	# Store the longs and lats of the antennas in lists - convert to
	# canonical units.

	longs = {}
	lats = {}
	radii = {}

	for ant in name:

		position = antenna_position[ant]
		radius = position['m2']['value']
		radius_unit = position['m2']['unit']
		radius_quantum = qaLoc.quantity( radius, radius_unit )
		radius_quantum = qaLoc.convert( radius_quantum, 'm' )
		radius = qaLoc.getvalue( radius_quantum )

		long = position['m0']['value']
		long_unit = position['m0']['unit']
		long_quantum = qaLoc.quantity( long, long_unit )
		long_quantum = qaLoc.convert( long_quantum, 'rad' )
		long = qaLoc.getvalue( long_quantum )

		lat = position['m1']['value']
		lat_unit = position['m1']['unit']
		lat_quantum = qaLoc.quantity( lat, lat_unit )
		lat_quantum = qaLoc.convert( lat_quantum, 'rad' )
		lat = qaLoc.getvalue( lat_quantum )

		longs[ant] = long
		lats[ant] = lat
		radii[ant] = radius


	# Get a median 'centre' for the array and derive x, y relative
	# to that.

	long_vals = np.array( longs.values() )
	long_vals -= np.median( long_vals )

	lat_vals = np.array( lats.values() )
	radii_vals = np.array( radii.values() )

	names = longs.keys()


	# Multiply longs by cos(lat) and radius to convert to metres.

	x = long_vals * np.cos(lat_vals) * radii_vals
	y = lat_vals * radii_vals


	# Make x,y relative to 'centre' of array.

	x -= np.median(x)
	y -= np.median(y)

	antenna_distance = {}

	for i,ant in enumerate(names):
		antenna_distance[ant] = np.sqrt( pow(x[i],2) + pow(y[i],2) )

	antennaInfoRead = True


        # Assign points to antennas based on their distance from the array
        # centre.

	distance_vals = antenna_distance.values()
	names = antenna_distance.keys()
	distance_argsort = np.argsort( np.array(distance_vals) )
	distance_points = {}

	for k,i in enumerate( distance_argsort ):
		ant = names[i]
		distance_points[ant] = len(distance_argsort) - k


	# Get the number of unflagged rows for each antenna in these spectral
	# windows / fields.

	tbLoc.open( '%s' % ms )

        taql = '(FIELD_ID IN %s)' % fieldList
	taql += ' && (DATA_DESC_ID IN %s)' % spwList
	taql += ' && (NOT(FLAG_ROW))'

	nrows = {}
	ids = np.arange( len(name) )

	for ant in name:
		ant_id = ids[name==ant][0]
		subTable = tbLoc.query(
		  '%s && (ANTENNA1==%s || ANTENNA2==%s)' % (taql,ant_id,ant_id)
		)
		nrows[ant] = subTable.nrows()
		subTable.close()

        tbLoc.close()


        # Assign points to antennas for the amount of data they have.

        nrows_vals = np.array( nrows.values(), np.float )
        names = np.array( nrows.keys() )
        nrows_normalised = (nrows_vals / max(nrows_vals)) * len(nrows_vals)
        nrows_points_array = map( int, nrows_normalised )

        nrows_points = {}

        for i,ant in enumerate(names):
		nrows_points[ant] = nrows_points_array[i]


        # Add up points for antennas.

	total_points = {}

	for ant in name:
		total_points[ant] = distance_points[ant] + nrows_points[ant]


        # Select antenna with the highest total score.

	total_points_array = np.array( total_points.values() )
	total_points_array = - total_points_array
	merit_argsort = np.argsort( total_points_array )

	sorted_antennas = names[merit_argsort]
	reference_antenna = names[merit_argsort[0]]


	# Return the reference antenna, sorted antenna list, and points

	return reference_antenna, sorted_antennas, total_points, \
	    distance_points, nrows_points

######################################################################

def uniq(inlist):
   uniques = []
   for item in inlist:
      if item not in uniques:
         uniques.append(item)
   return uniques

######################################################################

def find_EVLA_band(frequency):
   FLOW = [ 0.0e6, 150.0e6, 700.0e6, 2.0e9, 4.0e9, 8.0e9, 12.0e9, 18.0e9, 26.5e9, 40.0e9 ]
   FHIGH = [ 150.0e6, 700.0e6, 2.0e9, 4.0e9, 8.0e9, 12.0e9, 18.0e9, 26.5e9, 40.0e9, 56.0e9 ]
   BBAND = [ '4', 'P', 'L', 'S', 'C', 'X', 'U', 'K', 'A', 'Q' ]

   band = '?'
   for ii in range(0,len(FLOW)):
      if ((frequency > FLOW[ii]) and (frequency <= FHIGH[ii])):
         band = BBAND[ii]
   return band

######################################################################

def find_standards(positions):
# set the max separation as ~1'
   MAX_SEPARATION = 60*2.0e-5
   position_3C48 = me.direction('j2000', '1h37m41.299', '33d9m35.133')
   fields_3C48 = []
   position_3C138 = me.direction('j2000', '5h21m9.886', '16d38m22.051')
   fields_3C138 = []
   position_3C147 = me.direction('j2000', '5h42m36.138', '49d51m7.234')
   fields_3C147 = []
   position_3C286 = me.direction('j2000', '13h31m8.288', '30d30m23.959')
   fields_3C286 = []

   for ii in range(0,len(positions)):
      position = me.direction('j2000', str(positions[ii][0])+'rad', str(positions[ii][1])+'rad')
      separation = me.separation(position,position_3C48)['value'] * pi/180.0
      if (separation < MAX_SEPARATION):
         fields_3C48.append(ii)
      else:
         separation = me.separation(position,position_3C138)['value'] * pi/180.0
         if (separation < MAX_SEPARATION):
            fields_3C138.append(ii)
         else:
            separation = me.separation(position,position_3C147)['value'] * pi/180.0
            if (separation < MAX_SEPARATION):
               fields_3C147.append(ii)
            else:
               separation = me.separation(position,position_3C286)['value'] * pi/180.0
               if (separation < MAX_SEPARATION):
                  fields_3C286.append(ii)

   fields = [ fields_3C48, fields_3C138, fields_3C147, fields_3C286 ]
   return fields

######################################################################

def getFlaggedSoln(calTable):
    '''
    This method will look at the specified calibration table and return the
    fraction of flagged solutions for each Antenna, SPW, Poln.  This assumes
    that the specified cal table will not have any channel dependent flagging.

    return structure is a dictionary with AntennaID and Spectral Window ID
    as the keys and returns a list of fractional flagging per polarization in
    the order specified in the Cal Table.

    Example:
    execfile('getFlaggedSoln.py')
    result = getFlaggedSoln('testdelay.k')
    result[17][4] = [1.0, 1.0]
    '''

    from taskinit import tbtool
    mytb = tbtool.create()

    mytb.open(calTable)
    antCol = mytb.getcol('ANTENNA1')
    flagCol = mytb.getcol('FLAG')
    calDescCol = mytb.getcol('CAL_DESC_ID')
    mytb.close()

    mytb.open(calTable+"/CAL_DESC")
    spwCol = mytb.getcol('SPECTRAL_WINDOW_ID')[0]
    mytb.done()
    
    # Initialize a list to hold the results
    # The order here is row[antenna][calDescCol]([polarization] for flags)
    rows = [[0] * (max(calDescCol)+1) for idx in range(max(antCol)+1)]
    flags = [[[0] * len(flagCol) for idx in range(max(calDescCol)+1)]
             for idx2 in range(max(antCol)+1)]

    # Ok now go through and for each antenna and calDesc get the number
    # of flagged and total number of entries
    # Here we assume that there is no channel dependent information
    for idx in range(len(antCol)):
        rows[antCol[idx]][calDescCol[idx]] += 1
        for poln in range(len(flagCol)):
            if flagCol[poln][0][idx]:
                flags[antCol[idx]][calDescCol[idx]][poln] += 1
                
    # Now create the output dictionary
    outDict = {}
    for antIdx in range(max(antCol)+1):
        outDict[antIdx] = {}
        for calDescIdx in range(max(calDescCol)+1):
            if rows[antIdx][calDescIdx] > 0:
                outDict[antIdx][spwCol[calDescIdx]] = \
                       [pF/float(rows[antIdx][calDescIdx]) \
                        for pF in flags[antIdx][calDescIdx]]



    return outDict
    
######################################################################

import math
def getOptimumSize(size):
    '''
    This method takes as input the a size parameter.  The return is the smallest
    integer Y which satisfies the following conditions:
    * Y > size
    * Y = 2^a*3^b*5^c where a,b, and c are non-negative integers and at least one
    of a or b is 0 and c is nonzero
    '''
    def evaluate(pow2, pow3, pow5):
        # Convience method to calculate the value given multiples
        return int(math.pow(2,pow2) *math.pow(3,pow3)*math.pow(5,pow5))
    
    max5 = int(math.ceil(math.log(size,5)))
    returnValue = evaluate(0, 0, max5)
    for pow5 in range(max5,0,-1):
        pow2 = math.ceil(math.log(size/math.pow(5,pow5),2))
        if not pow2 < 0:
            returnValue = min(returnValue, evaluate(pow2,0,pow5))

        pow3 = math.ceil(math.log(size/math.pow(5,pow5),3))
        if not pow3 < 0:
            returnValue = min(returnValue, evaluate(0,pow3,pow5))
    return returnValue

######################################################################

import urllib2
import datetime

def correct_ant_posns (vis_name, print_offsets=False):
    '''
    Given an input visibility MS name (vis_name), find the antenna
    position offsets that should be applied.  This application should
    be via the gencal task, using caltype='antpos'.

    If the print_offsets parameter is True, will print out each of
    the found offsets (or indicate that none were found), otherwise
    runs silently.

    A list is returned where the first element is the returned error
    code, the second element is a string of the antennas, and the 
    third element is a list of antenna Bx,By,Bz offsets.  An example 
    return list might look like:
    [ 0, 'ea01,ea19', [0.0184, -0.0065, 0.005, 0.0365, -0.0435, 0.0543] ]

    Usage examples:

       CASA <1>: antenna_offsets = correct_ant_posns('test.ms')
       CASA <2>: if (antenna_offsets[0] == 0):
       CASA <3>:     gencal(vis='test.ms', caltable='cal.G', \
                     caltype='antpos', antenna=antenna_offsets[1], \
                     parameter=antenna_offsets[2])

    This function does NOT work for VLA datasets, only EVLA.  If an
    attempt is made to use the function for VLA data (prior to 2010),
    an error code of 1 is returned.

    The offsets are retrieved over the internet.  A description and the
    ability to manually examine and retrieve offsets is at:
    http://www.vla.nrao.edu/astro/archive/baselines/
    If the attempt to establish the internet connection fails, an error
    code of 2 is returned.

    Uses the same algorithm that the AIPS task VLANT does.


    bjb
    nrao
    spring 2012
    '''

    MONTHS = [ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' ]
    URL_BASE = 'http://www.vla.nrao.edu/cgi-bin/evlais_blines.cgi?Year='

#
# get start date+time of observation
#
    observation = tb.open(vis_name+'/OBSERVATION')
    time_range = tb.getcol('TIME_RANGE')
    tb.close()
    MJD_start_time = time_range[0][0] / 86400
    q1 = qa.quantity(time_range[0][0],'s')
    date_time = qa.time(q1,form='ymd')
# date_time looks like: '2011/08/10/06:56:49'
    [obs_year,obs_month,obs_day,obs_time_string] = date_time[0].split('/')
    if (int(obs_year) < 2010):
        if (print_offsets):
            print 'Does not work for VLA observations'
        return [1, '', []]
    [obs_hour,obs_minute,obs_second] = obs_time_string.split(':')
    obs_time = 10000*int(obs_year) + 100*int(obs_month) + int(obs_day) + \
               int(obs_hour)/24.0 + int(obs_minute)/1440.0 + \
               int(obs_second)/86400.0

#
# get antenna to station mappings
#
    observation = tb.open(vis_name+'/ANTENNA')
    ant_names = tb.getcol('NAME')
    ant_stations = tb.getcol('STATION')
    tb.close()
    ant_num_stas = []
    for ii in range(len(ant_names)):
        ant_num_stas.append([int(ant_names[ii][2:]), ant_names[ii], \
                            ant_stations[ii], 0.0, 0.0, 0.0, False])

    correction_lines = []
    current_year = datetime.datetime.now().year
# first, see if the internet connection is possible
    try:
        response = urllib2.urlopen(URL_BASE + '2010')
    except URLError, err:
        if (print_offsets):
            print 'No internet connection to antenna position correction URL ', \
                  err.reason
        return [2, '', []]
    response.close()
    for year in range(2010,current_year+1):
        response = urllib2.urlopen(URL_BASE + str(year))
        html = response.read()
        response.close()
        html_lines = html.split('\n')

        for correction_line in html_lines:
            if len(correction_line) and correction_line[0] != '<' and correction_line[0] != ';':
                for month in MONTHS:
                    if month in correction_line:
                        correction_lines.append(str(year)+' '+correction_line)
                        break

    corrections_list = []
    for correction_line in correction_lines:
        correction_line_fields = correction_line.split()
        if (len(correction_line_fields) > 9):
            [c_year, moved_date, obs_date, put_date, put_time_str, ant, pad, Bx, By, Bz] = correction_line_fields
            s_moved = moved_date[:3]
            i_month = 1
            for month in MONTHS:
                if (moved_date.find(month) >= 0):
                    break
                i_month = i_month + 1
            moved_time = 10000 * int(c_year) + 100 * i_month + \
                         int(moved_date[3:])
        else:
            [c_year, obs_date, put_date, put_time_str, ant, pad, Bx, By, Bz] = correction_line_fields
            moved_date = '     '
            moved_time = 0
        s_obs = obs_date[:3]
        i_month = 1
        for month in MONTHS:
            if (s_obs.find(month) >= 0):
                break
            i_month = i_month + 1
        obs_time_2 = 10000 * int(c_year) + 100 * i_month + int(obs_date[3:])
        s_put = put_date[:3]
        i_month = 1
        for month in MONTHS:
            if (s_put.find(month) >= 0):
                break
            i_month = i_month + 1
        put_time = 10000 * int(c_year) + 100 * i_month + int(put_date[3:])
        [put_hr, put_min] = put_time_str.split(':')
        put_time += (int(put_hr)/24.0 + int(put_min)/1440.0)
        corrections_list.append([c_year, moved_date, moved_time, obs_date, obs_time_2, put_date, put_time, int(ant), pad, float(Bx), float(By), float(Bz)])

    for correction_list in corrections_list:
        [c_year, moved_date, moved_time, obs_date, obs_time_2, put_date, put_time, ant, pad, Bx, By, Bz] = correction_list
        ant_ind = -1
        for ii in range(len(ant_num_stas)):
            ant_num_sta = ant_num_stas[ii]
            if (ant == ant_num_sta[0]):
                ant_ind = ii
                break
        if ((ant_ind == -1) or (ant_num_sta[6])):
# the antenna in this correction isn't in the observation, or is done, 
# so skip it
            pass
        ant_num_sta = ant_num_stas[ant_ind]
        if (moved_time):
# the antenna moved
            if (moved_time > obs_time):
# we are done considering this antenna
                ant_num_sta[6] = True
            else:
# otherwise, it moved, so the offsets should be reset
                ant_num_sta[3] = 0.0
                ant_num_sta[4] = 0.0
                ant_num_sta[5] = 0.0
        if ((put_time > obs_time) and (not ant_num_sta[6]) and \
            (pad == ant_num_sta[2])):
# it's the right antenna/pad; add the offsets to those already accumulated
            ant_num_sta[3] += Bx
            ant_num_sta[4] += By
            ant_num_sta[5] += Bz

    ants = []
    parms = []
    for ii in range(len(ant_num_stas)):
        ant_num_sta = ant_num_stas[ii]
        if ((ant_num_sta[3] != 0.0) or (ant_num_sta[4] != 0.0) or \
            (ant_num_sta[3] != 0.0)):
            if (print_offsets):
                print "offsets for antenna %4s : %8.5f  %8.5f  %8.5f" % \
                      (ant_num_sta[1], ant_num_sta[3], ant_num_sta[4], ant_num_sta[5])
            ants.append(ant_num_sta[1])
            parms.append(ant_num_sta[3])
            parms.append(ant_num_sta[4])
            parms.append(ant_num_sta[5])
    if ((len(parms) == 0) and print_offsets):
        print "No offsets found for this MS"
    ant_string = ','.join(["%s" % ii for ii in ants])
    return [ 0, ant_string, parms ]

######################################################################

def spwsforfield(vis,field):
    # get observed DDIDs for specified field from MAIN
    tb.open(vis)
    st=tb.query('FIELD_ID=='+str(field))
    ddids=pl.unique(st.getcol('DATA_DESC_ID'))
    st.close()
    tb.close()

    # get SPW_IDs corresponding to those DDIDs
    tb.open(vis+'/DATA_DESCRIPTION')
    spws=tb.getcol('SPECTRAL_WINDOW_ID')[ddids]
    tb.close()
    
    # return as a list
    return list(spws)

######################################################################

# ------------------------------------------------------------------------------

# findrefant.py

# Description:
# ------------
# This file contains the reference antenna heuristics.

# The present heuristics are geometry and flagging.

# Classes:
# --------
# RefAntHeuristics - This class chooses the reference antenna heuristics.
# RefAntGeometry   - This class contains the geometry heuristics for the
#                    reference antenna.
# RefAntFlagging   - This class contains the flagging heuristics for the
#                    reference antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.
# 2012 Jun 06 - Nick Elias, NRAO
#               Modified to exclude ALMA pipeline classe dependencies.

# ------------------------------------------------------------------------------

# Imports
# -------

import numpy

# ------------------------------------------------------------------------------
# class RefAntHeuristics
# ------------------------------------------------------------------------------

# RefAntHeuristics
# ----------------

# Description:
# ------------
# This class chooses the reference antenna heuristics.

# Public member variables:
# ------------------------
# vis      - This python string contains the MS name.
#
# field    - This python string or list of strings contains the field numbers
#            or IDs.  Presently it is used only for the flagging heuristic.
# spw      - This python string or list of strings contains the spectral
#            window numbers of IDs.  Presently it is used only for the
#            flagging heuristic.
# intent   - This python string or list of strings contains the intent(s).
#            Presently it is used only for the flagging heuristic.
#
# geometry - This python boolean determines whether the geometry heuristic will
#            be used.
# flagging - This python boolean determines whether the flagging heuristic will
#            be used.

# Public member functions:
# ------------------------
# __init__  - This public member function constructs an instance of the
#             RefAntHeuristics() class.
# calculate - This public member function forms the reference antenna list
#             calculated from the selected heuristics.

# Private member functions:
# -------------------------
# _get_names - This private member function gets the antenna names from the MS.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version created with public member variables vis, field,
#               spw, intent, geometry, and flagging; public member functions
#               __init__() and calculate(); and private member function
#               _get_names().
# 2012 Jun 06 - Nick Elias, NRAO
#               api inheritance eliminated.

# ------------------------------------------------------------------------------

class RefAntHeuristics:

# ------------------------------------------------------------------------------

# RefAntHeuristics::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntHeuristics()
# class.

# NB: If all of the defaults are chosen, no reference antenna list is returned.

# Inputs:
# -------
# vis        - This python string contains the MS name.
#
# field      - This python string or list of strings contains the field numbers
#              or IDs.  Presently it is used only for the flagging heuristic.
#              The default is ''.
# spw        - This python string or list of strings contains the spectral
#              window numbers of IDs.  Presently it is used only for the
#              flagging heuristic.  The default is ''.
# intent     - This python string or list of strings contains the intent(s).
#              Presently it is used only for the flagging heuristic.  The
#              default is ''.
#
# geometry   - This python boolean determines whether the geometry heuristic
#              will be used in automatic mode.  The default is False.
# flagging   - This python boolean determines whether the flagging heuristic
#              will be used in automatic mode.  The default is False.

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.
# 2012 Jun 06 - Nick Eluas, NRAO
#               Input parameter defaults added.

# ------------------------------------------------------------------------------

	def __init__( self, vis, field='', spw='', intent='', geometry=False,
	    flagging=False ):

		# Initialize the public member variables of this class

		self.vis = vis

		self.field = field
		self.spw = spw
		self.intent = intent

		self.geometry = geometry
		self.flagging = flagging


		# Return None

		return None

# ------------------------------------------------------------------------------

# RefAntHeuristics::calculate

# Description:
# ------------
# This public member function forms the reference antenna list calculated from
# the selected heuristics.

# NB: A total score is calculated from all heuristics.  The best antennas have
# the highest scores, so a reverse sort is performed to obtain the final list.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The numpy array of strings containing the ranked reference antenna list,
# returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def calculate( self ):

		# If no heuristics are specified, return no reference antennas

		if not ( self.geometry or self.flagging ): return []


		# Get the antenna names and initialize the score dictionary

		names = self._get_names()

		score = dict()
		for n in names: score[n] = 0.0


		# For each selected heuristic, add the score for each antenna

		self.geoScore = 0.0
		self.flagScore = 0.0

		if self.geometry:
			geoClass = RefAntGeometry( self.vis )
			self.geoScore = geoClass.calc_score()
			for n in names: score[n] += self.geoScore[n]

		if self.flagging:
			flagClass = RefAntFlagging( self.vis, self.field,
			    self.spw, self.intent )
			self.flagScore = flagClass.calc_score()
			for n in names:
                            try:
                                score[n] += self.flagScore[n]
                            except KeyError, e:
                                logprint ("WARNING: antenna "+str(e)+", is completely flagged and missing from calibrators.ms", logfileout='logs/refantwarnings.log')


		# Calculate the final score and return the list of ranked
		# reference antennas.  NB: The best antennas have the highest
		# score, so a reverse sort is required.

		keys = numpy.array( score.keys() )
		values = numpy.array( score.values() )
		argSort = numpy.argsort( values )[::-1]

		refAntUpper = keys[argSort]

		refAnt = list()
		for r in refAntUpper: refAnt.append( r.lower() )


		# Return the list of ranked reference antennas

		return( refAnt )

# ------------------------------------------------------------------------------

# RefAntHeuristics::_get_names

# Description:
# ------------
# This private member function gets the antenna names from the MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The numpy array of strings containing the antenna names, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_names( self ):

		# Create the local instance of the table tool and open the MS

		tbLoc = casac.table()
		tbLoc.open( self.vis+'/ANTENNA' )


		# Get the antenna names and capitalize them (unfortunately,
		# some CASA tools capitalize them and others don't)

		names = tbLoc.getcol( 'NAME' ).tolist()

		rNames = range( len(names) )
		for n in rNames: names[n] = names[n]


		# Close the local instance of the table tool and delete it

		tbLoc.close()
		del tbLoc


		# Return the antenna names

		return names

# ------------------------------------------------------------------------------
# class RefAntGeometry
# ------------------------------------------------------------------------------

# RefAntGeometry
# --------------

# Description:
# ------------
# This class contains the geometry heuristics for the reference antenna.

# Algorithm:
# ----------
# * Calculate the antenna distances from the array center.
# * Normalize the distances by the maximum distance.
# * Calculate the score for each antenna, which is one minus the normalized
#   distance.  The best antennas have the highest score.
# * Sort according to score.

# Public member variables:
# ------------------------
# vis - This python string contains the MS name.

# Public member functions:
# ------------------------
# __init__   - This public member function constructs an instance of the
#              RefAntGeometry() class.
# calc_score - This public member function calculates the geometry score for
#              each antenna.

# Private member functions:
# -------------------------
# _get_info       - This private member function gets the information from the
#                   antenna table of the MS.
# _get_measures   - This private member function gets the measures from the
#                   antenna table of the MS.
# _get_latlongrad - This private member function gets the latitude, longitude
#                   and radius (from the center of the earth) for each antenna.
# _calc_distance  - This private member function calculates the antenna
#                   distances from the array reference from the radii,
#                   longitudes, and latitudes.
# _calc_score     - This private member function calculates the geometry score
#                   for each antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

class RefAntGeometry:

# ------------------------------------------------------------------------------

# RefAntGeometry::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntGeometry()
# class.

# Inputs:
# -------
# vis - This python string contains the MS name.

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def __init__( self, vis ):

		# Set the public variables

		self.vis = vis


		# Return None

		return None

# ------------------------------------------------------------------------------

# RefAntGeometry::calc_score

# Description:
# ------------
# This public member function calculates the geometry score for each antenna.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def calc_score( self ):

		# Get the antenna information, measures, and locations

		info = self._get_info()
		measures = self._get_measures( info )
		radii, longs, lats = self._get_latlongrad( info, measures )


		# Calculate the antenna distances and scores

		distance = self._calc_distance( radii, longs, lats )
		score = self._calc_score( distance )


		# Return the scores

		return score

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_info

# Description:
# ------------
# This private member function gets the information from the antenna table of
# the MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the antenna information, returned via the
# function value.  The dictionary format is:
# 'position'          - This numpy array contains the antenna positions.
# 'flag_row'          - This numpy array of booleans contains the flag row
#                       booleans.  NB: This element is of limited use now and
#                       may be eliminated.
# 'name'              - This numpy array of strings contains the antenna names.
# 'position_keywords' - This python dictionary contains the antenna information.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_info( self ):

		# Create the local instance of the table tool and open it with
		# the antenna subtable of the MS

		tbLoc = casac.table()
		tbLoc.open( self.vis+'/ANTENNA' )


		# Get the antenna information from the antenna table

		info = dict()

		info['position'] = tbLoc.getcol( 'POSITION' )
		info['flag_row'] = tbLoc.getcol( 'FLAG_ROW' )
		info['name'] = tbLoc.getcol( 'NAME' )
		info['position_keywords'] = tbLoc.getcolkeywords( 'POSITION' )


		# Close the table tool and delete the local instance

		tbLoc.close()
		del tbLoc


		# The flag tool appears to return antenna names as upper case,
		# which seems to be different from the antenna names stored in
		# MSes.  Therefore, these names will be capitalized here.

		rRow = range( len( info['name'] ) )


		# Return the antenna information

		return info

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_measures

# Description:
# ------------
# This private member function gets the measures from the antenna table of the
# MS.

# Inputs:
# -------
# info - This python dictionary contains the antenna information from private
#        member function _get_info().

# Outputs:
# --------
# The python dictionary containing the antenna measures, returned via the
# function value.  The dictionary format is:
# '<antenna name>' - The python dictionary containing the antenna measures.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_measures( self, info ):

		# Create the local instances of the measures and quanta tools

		meLoc = casac.measures()
		qaLoc = casac.quanta()


		# Initialize the measures dictionary and the position and
		# position_keywords variables

		measures = dict()

		position = info['position']
		position_keywords = info['position_keywords']

		rf = position_keywords['MEASINFO']['Ref']

		for row,ant in enumerate( info['name'] ):

			if not info['flag_row'][row]:

				p = position[0,row]
				pk = position_keywords['QuantumUnits'][0]
				v0 = qaLoc.quantity( p, pk )

				p = position[1,row]
				pk = position_keywords['QuantumUnits'][1]
				v1 = qaLoc.quantity( p, pk )

				p = position[2,row]
				pk = position_keywords['QuantumUnits'][2]
				v2 = qaLoc.quantity( p, pk )

				measures[ant] = meLoc.position( rf=rf, v0=v0,
				    v1=v1, v2=v2 )


		# Delete the local instances of the measures and quanta tools

		del qaLoc
		del meLoc


		# Return the measures

		return measures

# ------------------------------------------------------------------------------

# RefAntGeometry::_get_latlongrad

# Description:
# ------------
# This private member function gets the latitude, longitude and radius (from the
# center of the earth) for each antenna.

# Inputs:
# -------
# info     - This python dictionary contains the antenna information from
#            private member function _get_info().
# measures - This python dictionary contains the antenna measures from private
#            member function _get_measures().

# Outputs:
# --------
# The python tuple containing containing radius, longitude, and latitude python
# dictionaries, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_latlongrad( self, info, measures ):

		# Create the local instance of the quanta tool

		qaLoc = casac.quanta()


		# Get the radii, longitudes, and latitudes

		radii = dict()
		longs = dict()
		lats = dict()

		for ant in info['name']:

			value = measures[ant]['m2']['value']
			unit = measures[ant]['m2']['unit']
			quantity = qaLoc.quantity( value, unit )
			convert = qaLoc.convert( quantity, 'm' )
			radii[ant] = qaLoc.getvalue( convert )

			value = measures[ant]['m0']['value']
			unit = measures[ant]['m0']['unit']
			quantity = qaLoc.quantity( value, unit )
			convert = qaLoc.convert( quantity, 'rad' )
			longs[ant] = qaLoc.getvalue( convert )

			value = measures[ant]['m1']['value']
			unit = measures[ant]['m1']['unit']
			quantity = qaLoc.quantity( value, unit )
			convert = qaLoc.convert( quantity, 'rad' )
			lats[ant] = qaLoc.getvalue( convert )


		# Delete the local instance of the quanta tool

		del qaLoc


		# Return the tuple containing the radius, longitude, and
		# latitude python dictionaries

		return radii, longs, lats

# ------------------------------------------------------------------------------

# RefAntGeometry::_calc_distance

# Description:
# ------------
# This private member function calculates the antenna distances from the array
# reference from the radii, longitudes, and latitudes.

# NB: The array reference is the median location.

# Inputs:
# -------
# radii - This python dictionary contains the radius (from the center of the
#         earth) for each antenna.
# longs - This python dictionary contains the longitude for each antenna.
# lats  - This python dictionary contains the latitude for each antenna.

# Outputs:
# --------
# The python dictionary containing the antenna distances from the array
# reference, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _calc_distance( self, radii, longs, lats ):

		# Convert the dictionaries to numpy float arrays.  The median
		# longitude is subtracted.

		radiusValues = numpy.array( radii.values() )

		longValues = numpy.array( longs.values() )
		longValues -= numpy.median( longValues )

		latValues = numpy.array( lats.values() )


		# Calculate the x and y antenna locations.  The medians are
		# subtracted.

		x = longValues * numpy.cos(latValues) * radiusValues
		x -= numpy.median( x )

		y = latValues * radiusValues
		y -= numpy.median( y )


		# Calculate the antenna distances from the array reference and
		# return them

		distance = dict()
		names = radii.keys()

		for i,ant in enumerate(names):
			distance[ant] = numpy.sqrt( pow(x[i],2) + pow(y[i],2) )

		return distance

# ------------------------------------------------------------------------------

# RefAntGeometry::_calc_score

# Description:
# ------------
# This private member function calculates the geometry score for each antenna.

# Algorithm:
# ----------
# * Calculate the antenna distances from the array center.
# * Normalize the distances by the maximum distance.
# * Calculate the score for each antenna, which is one minus the normalized
#   distance.  The best antennas have the highest score.
# * Sort according to score.

# Inputs:
# -------
# distance - This python dictionary contains the antenna distances from the
#            array reference.  They are calculated in private member function
#            _calc_distance().

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _calc_score( self, distance ):

		# Get the number of good data, calculate the fraction of good
		# data, and calculate the good and bad weights

		far = numpy.array( distance.values(), numpy.float )
		fFar = far / float( numpy.max(far) )

		wFar = fFar * len(far)
		wClose = ( 1.0 - fFar ) * len(far)


		# Calculate the score for each antenna and return them

		score = dict()

		names = distance.keys()
		rName = range( len(wClose) )

		for n in rName: score[names[n]] = wClose[n]

		return score

# ------------------------------------------------------------------------------

# RefAntFlagging
# --------------

# Description:
# ------------
# This class contains the flagging heuristics for the reference antenna.

# Algorithm:
# ----------
# * Get the number of unflagged (good) data for each antenna.
# * Normalize the good data by the maximum good data.
# * Calculate the score for each antenna, which is one minus the normalized
#   number of good data.  The best antennas have the highest score.
# * Sort according to score.

# Public member variables:
# ------------------------
# vis    - This python string contains the MS name.
#
# field  - This python string or list of strings contains the field numbers or
#          or IDs.
# spw    - This python string or list of strings contains the spectral window
#          numbers of IDs.
# intent - This python string or list of strings contains the intent(s).

# Public member functions:
# ------------------------
# __init__   - This public member function constructs an instance of the
#              RefAntFlagging() class.
# calc_score - This public member function calculates the flagging score for
#              each antenna.

# Private member functions:
# -------------------------
# _get_good   - This private member function gets the number of unflagged (good)
#               data from the MS.
# _calc_score - This private member function calculates the flagging score for
#               each antenna.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

class RefAntFlagging:

# ------------------------------------------------------------------------------

# RefAntFlagging::__init__

# Description:
# ------------
# This public member function constructs an instance of the RefAntFlagging()
# class.

# Inputs:
# -------
# vis    - This python string contains the MS name.
#
# field  - This python string or list of strings contains the field numbers or
#          or IDs.
# spw    - This python string or list of strings contains the spectral window
#          numbers of IDs.
# intent - This python string or list of strings contains the intent(s).

# Outputs:
# --------
# None, returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def __init__( self, vis, field, spw, intent ):

		# Set the public member functions

		self.vis = vis

		self.field = field
		self.spw = spw
		self.intent = intent


		# Return None

		return None

# ------------------------------------------------------------------------------

# RefAntFlagging::calc_score

# Description:
# ------------
# This public member function calculates the flagging score for each antenna.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def calc_score( self ):

		# Calculate the number of unflagged (good) measurements for each
		# antenna, determine the score, and return them

		good = self._get_good()
		score = self._calc_score( good )

		return( score )

# ------------------------------------------------------------------------------

# RefAntFlagging::_get_good

# Description:
# ------------
# This private member function gets the number of unflagged (good) data from the
# MS.

# Inputs:
# -------
# None.

# Outputs:
# --------
# The dictionary containing the number of unflagged (good) data from the MS,
# returned via the function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _get_good( self ):

                #Update April 2015 to use the flagging task
                task_args = {'vis'          : self.vis, 
                           'mode'         : 'summary',
                           'field'        : self.field,
                           'spw'          : self.spw,
                           'intent'       : self.intent,
                           'display'      : '',
                           'flagbackup'   : False,
                           'savepars'     : False}

                d = flagdata(**task_args)                                               



		# Calculate the number of good data for each antenna and return
		# them

		antenna = d['antenna']
		good = dict()

		for a in antenna.keys():
			good[a] = antenna[a]['total'] - antenna[a]['flagged']

		return( good )

# ------------------------------------------------------------------------------

# RefAntFlagging::_calc_score

# Description:
# ------------
# This private member function calculates the flagging score for each antenna.

# Algorithm:
# ----------
# * Get the number of unflagged (good) data for each antenna.
# * Normalize the good data by the maximum good data.
# * Calculate the score for each antenna, which is one minus the normalized
#   number of good data.  The best antennas have the highest score.
# * Sort according to score.

# Inputs:
# -------
# good - This python dictionary contains the number of unflagged (good) data
#        from the MS.  They are obtained in private member function _get_good().

# Outputs:
# --------
# The python dictionary containing the score for each antenna, returned via the
# function value.

# Modification history:
# ---------------------
# 2012 May 21 - Nick Elias, NRAO
#               Initial version.

# ------------------------------------------------------------------------------

	def _calc_score( self, good ):

		# Get the number of good data, calculate the fraction of good
		# data, and calculate the good and bad weights

		nGood = numpy.array( good.values(), numpy.float )
		fGood = nGood / float( numpy.max(nGood) )

		wGood = fGood * len(nGood)
		wBad = ( 1.0 - fGood ) * len(nGood)


		# Calculate the score for each antenna and return them

		score = dict()

		names = good.keys()
		rName = range( len(wGood) )

		for n in rName: score[names[n]] = wGood[n]

		return score

######################################################################

def testBPdgains(calMs,calTable,calSpw,calScans,calSolint,refAnt,minBL_for_cal,priorcals,do3C84,UVrange3C84):

    GainTables=copy.copy(priorcals)
    GainTables.append('testdelay.k')

    if (do3C84 == True):
        gaincal(vis=calMs,caltable=calTable,field='',spw=calSpw,intent='',\
            selectdata=True,uvrange=UVrange3C84,scan=calScans,\
            solint=calSolint,combine='scan',\
            preavg=-1.0,refant=refAnt,minblperant=minBL_for_cal,minsnr=5.0,\
            solnorm=False,\
            gaintype='G',smodel=[],calmode='ap',append=False,\
            docallib=False,\
            gaintable=GainTables,gainfield=[''],interp=[''],\
            spwmap=[],parang=False)
    else:
        gaincal(vis=calMs,caltable=calTable,field='',spw=calSpw,intent='',\
            selectdata=True,uvrange='',scan=calScans,\
            solint=calSolint,combine='scan',\
            preavg=-1.0,refant=refAnt,minblperant=minBL_for_cal,minsnr=5.0,\
            solnorm=False,\
            gaintype='G',smodel=[],calmode='ap',append=False,\
            docallib=False,\
            gaintable=GainTables,gainfield=[''],interp=[''],\
            spwmap=[],parang=False)

    flaggedSolnResult=getCalFlaggedSoln(calTable)

    return flaggedSolnResult

######################################################################

def testdelays(calMs,calTable,calField,calScans,refAnt,minBL_for_cal,priorcals,do3C84,UVrange3C84):

    GainTables=copy.copy(priorcals)
    GainTables.append('testdelayinitialgain.g')

# NB: can't use uvrange for delay cals because it flags all antennas
# beyond uvrange from the refant; so leave as all in delay cal and
# have any residual structure taken out by bandpass and gaincal
# until this is fixed

#    if (do3C84 == True):
#        gaincal(vis=calMs,caltable=calTable,field=calField,spw='',intent='',\
#            selectdata=True,uvrange=UVrange3C84,scan=calScans,solint='inf',\
#            combine='scan',preavg=-1.0,refant=refAnt,\
#            minblperant=minBL_for_cal,minsnr=3.0,solnorm=False,\
#            gaintype='K',smodel=[],calmode='p',append=False,\
#            docallib=False,\
#            gaintable=GainTables,gainfield=[''],\
#            interp=[''],spwmap=[],parang=False)
#    else:

    gaincal(vis=calMs,caltable=calTable,field=calField,spw='',intent='',\
        selectdata=True,uvrange='',scan=calScans,solint='inf',\
        combine='scan',preavg=-1.0,refant=refAnt,\
        minblperant=minBL_for_cal,minsnr=3.0,solnorm=False,\
        gaintype='K',smodel=[],calmode='p',append=False,\
        docallib=False,\
        gaintable=GainTables,gainfield=[''],\
        interp=[''],spwmap=[],parang=False)

    flaggedSolnResult=getCalFlaggedSoln(calTable)

    return flaggedSolnResult

######################################################################

def testgains(calMs,calTable,calSpw,calScans,calSolint,refAnt,minBL_for_cal,combtime):

# NB: when the pipeline can use the full MS instead of
# calibrators.ms, gaintable will have to include all priorcals plus
# delay.k and BPcal.b

    gaincal(vis=calMs,caltable=calTable,field='',spw=calSpw,intent='',\
        selectdata=True,scan=calScans,solint=calSolint,combine='scan',\
        preavg=-1.0,refant=refAnt,minblperant=minBL_for_cal,minsnr=5.0,\
        solnorm=False,\
        gaintype='G',smodel=[],calmode='ap',append=False,\
        docallib=False,\
        gaintable=[''],gainfield=[''],interp=[''],\
        spwmap=[],parang=False)

    flaggedSolnResult=getCalFlaggedSoln(calTable)

    return flaggedSolnResult

######################################################################

def semiFinaldelays(calMs,calTable,calField,calScans,refAnt,minBL_for_cal,priorcals,do3C84,UVrange3C84):

    GainTables=copy.copy(priorcals)
    GainTables.append('semiFinaldelayinitialgain.g')
    
# NB: can't use uvrange for delay cals because it flags all antennas
# beyond uvrange from the refant; so leave as all in delay cal and
# have any residual structure taken out by bandpass and gaincal
# until this is fixed

#    if (do3C84 == True):
#        gaincal(vis=calMs,caltable=calTable,field=calField,spw='',intent='',\
#            selectdata=True,uvrange=UVrange3C84,scan=calScans,solint='inf',\
#            combine='scan',preavg=-1.0,refant=refAnt,\
#            minblperant=minBL_for_cal,minsnr=3.0,solnorm=False,\
#            gaintype='K',smodel=[],calmode='p',append=False,\
#            docallib=False,\
#            gaintable=GainTables,gainfield=[''],\
#            interp=[''],spwmap=[],parang=False)
#    else:

    gaincal(vis=calMs,caltable=calTable,field=calField,spw='',intent='',\
        selectdata=True,uvrange='',scan=calScans,solint='inf',\
        combine='scan',preavg=-1.0,refant=refAnt,\
        minblperant=minBL_for_cal,minsnr=3.0,solnorm=False,\
        gaintype='K',smodel=[],calmode='p',append=False,\
        docallib=False,\
        gaintable=GainTables,gainfield=[''],\
        interp=[''],spwmap=[],parang=False)

    flaggedSolnResult=getCalFlaggedSoln(calTable)

    return flaggedSolnResult


######################################################################
#Added by B. Kent October 18, 2012

def find_3C84(positions):
    MAX_SEPARATION = 60*2.0e-5
    position_3C84 = me.direction('j2000', '3h19m48.160', '41d30m42.106')
    fields_3C84 = []

    for ii in range(0,len(positions)):
        position = me.direction('j2000', str(positions[ii][0])+'rad', str(positions[ii][1])+'rad')
        separation = me.separation(position,position_3C84)['value'] * pi/180.0
        if (separation < MAX_SEPARATION):
                  fields_3C84.append(ii)

    return fields_3C84

######################################################################

def checkblankplot(plotfile, maincasalog):
    ''' Returns True if file is removed
    '''
    blankplot = False
    fhandle = os.popen('tail ' + maincasalog + ' | grep "Plotting 0 unflagged points."')
    if fhandle.read() != '':
        #Plot is blank, so delete the file
        os.system('rm -rf ' + plotfile)
        logprint("Plot is blank - file removed.", logfileout='logs/targetflag.log')
        blankplot = True
#
    return blankplot

