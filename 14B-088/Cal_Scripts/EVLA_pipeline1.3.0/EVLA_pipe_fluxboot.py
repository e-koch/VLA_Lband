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

# DO THE FLUX DENSITY BOOTSTRAPPING

logprint("Starting EVLA_pipe_fluxboot.py", logfileout='logs/fluxboot.log')
time_list = runtiming('fluxboot', 'start')
QA2_fluxboot = 'Pass'

fluxscale_output = msname.rstrip('ms') + 'fluxdensities'
fluxcalfields = flux_field_select_string

# Check if there are other calibrators to bootstrap to. If not,
# we can skip this step.

if calibrator_field_select_string == flux_field_select_string:
    logprint("Flux density bootstrapping not needed. No other calibrators found",
             logfileout='logs/fluxboot.log')

else:
    logprint("Doing flux density bootstrapping",
             logfileout='logs/fluxboot.log')
    logprint("Flux densities will be written to " +
             fluxscale_output, logfileout='logs/fluxboot.log')

    syscommand = 'rm -rf ' + fluxscale_output
    os.system(syscommand)

    # If this is needed earlier in the script move to msinfo.py
    #currentcasalog = casalogger.func_globals['thelogfile']

    casalog.setlogfile(fluxscale_output)

    default('fluxscale')
    vis = 'calibrators.ms'
    caltable = 'fluxgaincal.g'
    fluxtable = 'fluxgaincalFcal.g'
    reference = [fluxcalfields]
    transfer = ['']
    listfile = ''
    append = False
    refspwmap = [-1]
    incremental = False
    fitorder = 1
    async = False
    fluxscale_result = fluxscale()

    casalog.setlogfile(maincasalog)

    logprint("Fitting data with power law", logfileout='logs/fluxboot.log')

    #
    # the variable center_frequencies should already have been filled out
    # with the reference frequencies of the spectral window table
    #

    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

    try:
        ff = open(fluxscale_output, 'r')
    except IOError as err:
        logprint(fluxscale_output + " doesn't exist, error: " +
                 err.filename, logfileout='logs/fluxboot.log')

    # looking for lines like:
    # 2012-03-09 21:30:23     INFO    fluxscale::::    Flux density for J1717-3342 in SpW=3 is: 1.94158 +/- 0.0123058 (SNR = 157.777, N= 34)
    # sometimes they look like:
    # 2012-03-09 21:30:23     INFO    fluxscale::::    Flux density for J1717-3342 in SpW=0 is:  INSUFFICIENT DATA
    # so watch for that.

    sources = []
    flux_densities = []
    spws = []

    # Find the field_ids in the dictionary returned from the CASA task
    # fluxscale
    dictkeys = fluxscale_result.keys()
    keys_to_remove = ['freq', 'spwName', 'spwID']
    dictkeys = [
        field_id for field_id in dictkeys if field_id not in keys_to_remove]

    for field_id in dictkeys:
        sourcename = fluxscale_result[field_id]['fieldName']
        secondary_keys = fluxscale_result[field_id].keys()
        secondary_keys_to_remove = [
            'fitRefFreq', 'spidxerr', 'spidx', 'fitFluxd', 'fieldName', 'fitFluxdErr']
        spwkeys = [
            spw_id for spw_id in secondary_keys if spw_id not in secondary_keys_to_remove]

        for spw_id in spwkeys:
            flux_d = list(fluxscale_result[field_id][spw_id]['fluxd'])
            flux_d_err = list(fluxscale_result[field_id][spw_id]['fluxdErr'])
            # spwslist  = list(int(spw_id))

            # flux_d = list(fluxscale_result[field_id]['fluxd'])
            # flux_d_err = list(fluxscale_result[field_id]['fluxdErr'])
            # spwslist  = list(fluxscale_result['spwID'])

            for i in range(0, len(flux_d)):
                if (flux_d[i] != -1.0 and flux_d[i] != 0.0):
                    sources.append(sourcename)
                    flux_densities.append(
                        [float(flux_d[i]), float(flux_d_err[i])])
                    spws.append(int(spw_id))

    ii = 0
    unique_sources = list(np.unique(sources))
    results = []
    for source in unique_sources:
        indices = []
        for ii in range(len(sources)):
            if (sources[ii] == source):
                indices.append(ii)
        bands = []
        for ii in range(len(indices)):
            bands.append(find_EVLA_band(center_frequencies[spws[indices[ii]]]))
        unique_bands = list(np.unique(bands))
        for band in unique_bands:
            lfreqs = []
            lfds = []
            lerrs = []
            uspws = []
            for ii in range(len(indices)):
                if find_EVLA_band(center_frequencies[spws[indices[ii]]]) == band:
                    lfreqs.append(log10(center_frequencies[spws[indices[ii]]]))
                    lfds.append(log10(flux_densities[indices[ii]][0]))
                    lerrs.append(
                        log10(e) * flux_densities[indices[ii]][1] / flux_densities[indices[ii]][0])
                    uspws.append(spws[indices[ii]])
    # if we didn't care about the errors on the data or the fit coefficients, just:
    #       coefficients = np.polyfit(lfreqs, lfds, 1)
    # or, if we ever get to numpy 1.7.x, for weighted fit, and returning
    # covariance matrix, do:
    #       ...
    #       weights = []
    #       weight_sum = 0.0
    #       for ii in range(len(lfreqs)):
    #           weights.append(1.0 / (lerrs[ii]*lerrs[ii]))
    #           weight_sum += weights[ii]
    #       for ii in range(len(weights)):
    #           weights[ii] /= weight_sum
    #       coefficients = np.polyfit(lfreqs, lfds, 1, w=weights, cov=True)
    # but, for now, use the full scipy.optimize.leastsq route...
    #
    # actually, after a lot of testing, np.polyfit does not return a global
    # minimum solution.  sticking with leastsq (modified as below to get the
    # proper errors), or once we get a modern enough version of scipy, moving
    # to curve_fit, is better.
    #

            if len(lfds) < 2:
                aa = lfds[0]
                bb = 0.0
                SNR = 0.0
            else:
                alfds = scp.array(lfds)
                alerrs = scp.array(lerrs)
                alfreqs = scp.array(lfreqs)
                pinit = [0.0, 0.0]
                fit_out = scpo.leastsq(
                    errfunc, pinit, args=(alfreqs, alfds, alerrs), full_output=1)
                pfinal = fit_out[0]
                covar = fit_out[1]
                aa = pfinal[0]
                bb = pfinal[1]
    #
    # the fit is of the form:
    #     log(S) = a + b * log(f)
    # with a = pfinal[0] and b = pfinal[1].  the errors on the coefficients are
    # sqrt(covar[i][i]*residual_variance) with the residual covariance calculated
    # as below (it's like the reduced chi squared without dividing out the errors).
    # see the scipy.optimize.leastsq documentation and
    # http://stackoverflow.com/questions/14854339/in-scipy-how-and-why-does-curve-fit-calculate-the-covariance-of-the-parameter-es
    #
                summed_error = 0.0
                for ii in range(len(alfds)):
                    model = aa + bb * alfreqs[ii]
                    residual = (model - alfds[ii]) * (model - alfds[ii])
                    summed_error += residual
                residual_variance = summed_error / (len(alfds) - 2)
                SNR = fabs(bb) / sqrt(covar[1][1] * residual_variance)

    #
    # take as the reference frequency the lowest one.  (this shouldn't matter,
    # in principle).
    #
            reffreq = 10.0**lfreqs[0] / 1.0e9
            fluxdensity = 10.0**(aa + bb * lfreqs[0])
            spix = bb
            results.append([source, uspws, fluxdensity, spix, SNR, reffreq])
            logprint(source + ' ' + band + ' fitted spectral index & SNR = ' +
                     str(spix) + ' ' + str(SNR), logfileout='logs/fluxboot.log')
            logprint("Frequency, data, error, and fitted data:",
                     logfileout='logs/fluxboot.log')
            for ii in range(len(lfreqs)):
                SS = fluxdensity * (10.0**lfreqs[ii] / reffreq / 1.0e9)**spix
                fderr = lerrs[ii] * (10**lfds[ii]) / log10(e)
                logprint('    ' + str(10.0**lfreqs[ii] / 1.0e9) + '  ' + str(10.0**lfds[
                         ii]) + '  ' + str(fderr) + '  ' + str(SS), logfileout='logs/fluxboot.log')

    logprint("Setting power-law fit in the model column",
             logfileout='logs/fluxboot.log')

    for result in results:
        for spw_i in result[1]:
    #
    # here, check on SNR, but don't do this yet, until we know what typical SNRs are
    #
    #       if result[4] > SNRlimit:
            logprint('Running setjy on spw ' + str(spw_i),
                     logfileout='logs/fluxboot.log')
            default('setjy')
            vis = 'calibrators.ms'
            field = str(result[0])
            #spw = ','.join(["%s" % ii for ii in result[1]])
            spw = str(spw_i)
            selectdata = False
            scalebychan = True
            standard = 'manual'
            fluxdensity = [result[2], 0, 0, 0]
            spix = result[3]
            reffreq = str(result[5]) + 'GHz'
            usescratch = False
            async = False
            setjy()
            vis = ms_active
            setjy()
            if (abs(spix) > 5.0):
                QA2_fluxboot = 'Fail'

    logprint("Flux density bootstrapping finished",
             logfileout='logs/fluxboot.log')

logprint("Plotting model calibrator flux densities",
         logfileout='logs/fluxboot.log')


syscommand = 'rm -rf bootstrappedFluxDensities.png'
os.system(syscommand)

# clearstat()

default('plotms')
vis = ms_active
xaxis = 'freq'
yaxis = 'amp'
ydatacolumn = 'model'
selectdata = True
scan = calibrator_scan_select_string
correlation = corrstring
averagedata = True
avgtime = '1e8s'
avgscan = True
transform = False
extendflag = False
iteraxis = ''
coloraxis = 'field'
plotrange = []
title = ''
xlabel = ''
ylabel = ''
showmajorgrid = False
showminorgrid = False
plotfile = 'bootstrappedFluxDensities.png'
overwrite = True
showgui = False
async = False
plotms()

logprint("QA2 score: " + QA2_fluxboot, logfileout='logs/fluxboot.log')
logprint("Finished EVLA_pipe_fluxboot.py", logfileout='logs/fluxboot.log')
time_list = runtiming('fluxboot', 'end')

pipeline_save()
