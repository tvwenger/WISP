"""
calibration.py - WISP Calibration Pipeline

Calibrate a measurement set by scripting CASA tasks and generating
diagnostic plots.

Copyright(C) 2018 by
Trey V. Wenger; tvwenger@gmail.com

GNU General Public License v3 (GNU GPLv3)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Changelog:
Trey V. Wenger November 2018 - V1.0

Trey V. Wenger December 2018 - V1.1
    Removed TFCROP - it was causing problems with VLA data.
    Added preliminary scan flagging, extend flags, and tunable
    parameters for shadow tolerance and quack interval.
    Added opacity, gaincurve, and antenna position calibrations.
    Changed plotcal figures to use plotms for generation.
"""

import __main__ as casa
import os
import numpy as np
import glob
import re
import time
import pickle
import logging
import logging.config
import ConfigParser
import gc

__version__ = "1.1"

# load logging configuration file
logging.config.fileConfig('logging.conf')

def natural_sort(l):
    """
    Natural sort an alphanumeric list

    Inputs: l
      l :: list of strings
        The list of strings to be sorted

    Returns: sorted_l
      sorted_l :: list of strings
        The sorted list
    """
    # Convert text to integers or lowercase
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    # define the sorting algorithm
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    # return the sorted list
    return sorted(l, key = alphanum_key)

def get_refant(vis=''):
    """
    Return the best reference antenna to use. For now, just prompt
    the user
    TODO: make this smarter but also efficient

    Inputs: string
      vis :: string
        The measurement set

    Returns: refant
      refant :: string
        The reference antenna
    """
    #casa.plotants(vis=vis)
    refant = raw_input("Refant? ")
    return refant

def setup(vis='',config=None):
    """
    Perform setup tasks: get reference antenna, generate listobs
    file, find line and continuum spectral windows, generate a list
    of primary calibrators, secondary calibrators, flux calibrators,
    and science targets, create needed directories.

    Inputs: vis, config
      vis :: string
        The measurement set
      config :: a ConfigParser objet
        The ConfigParser object for this project

    Returns: my_cont_spws, my_line_spws, flux_cals, primary_cals,
             secondary_cals, science_targets, refant
      my_cont_spws :: string
        comma-separated string of continuum spws
      my_line_spws :: string
        comma-separated string of line spws
      flux_cals :: list of strings
        list of flux calibrator names
      primary_cals :: list of strings
        list of primary calibrator names
      secondary_cals :: list of strings
        list of secondary calibrator names
      science_targets :: list of strings
        list of science targets names
      refant :: string
        The reference antenna
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None") 
    #
    # find good referance antenna
    #
    logger.info("Looking for good reference antenna...")
    refant = get_refant(vis=vis)
    if refant is None:
        logger.critical("Error: No good referance antenna found!")
        raise ValueError("No good referance antenna found!")
    logger.info("Done. Found reference antenna: {0}".format(refant))
    #
    # Generate listobs file
    #
    if not os.path.isfile('listobs.txt'):
        logger.info("Generating listobs file...")
        casa.listobs(vis=vis,listfile='listobs.txt')
        logger.info("Done.")
    #
    # Get continuum and line spws from configuration file
    #
    my_cont_spws = config.get("Spectral Windows", "Continuum")
    my_line_spws = config.get("Spectral Windows", "Line")
    logger.info("Found continuum spws: {0}".format(my_cont_spws))
    logger.info("Found line spws: {0}".format(my_line_spws))
    #
    # get field names
    #
    logger.info("Looking for field names...")
    fields = casa.vishead(vis=vis, mode='get', hdkey='field')[0]
    logger.info("Found fields:")
    logger.info('{0}'.format(fields))
    #
    # Get primary calibrator fields if they are not in config
    #
    if config.get('Calibrators', 'Primary Calibrators') == '':
        primary_cals = []
        logger.info("Looking for primary calibrators...")
        with open('listobs.txt') as f:
            for line in f:
                for field in fields:
                    if field in line and "CALIBRATE_BANDPASS" in line \
                      and field not in primary_cals:
                        primary_cals.append(field)
        logger.info("Done")
    else:
        primary_cals = [field for field in config.get('Calibrators', 'Primary Calibrators').splitlines()
                        if field in fields]
    logger.info("Primary calibrators: {0}".format(primary_cals))
    #
    # Get Secondary calibrator fields if they are not in config
    #
    if config.get('Calibrators', 'Secondary Calibrators') == '':
        secondary_cals = []
        logger.info("Looking for secondary calibrators...")
        with open('listobs.txt') as f:
            for line in f:
                for field in fields:
                    if field in line and ("CALIBRATE_AMPLI" in line or "CALIBRATE_PHASE" in line) \
                      and field not in secondary_cals:
                        secondary_cals.append(field)
        logger.info("Done")
    else:
        secondary_cals = [field for field in config.get('Calibrators','Secondary Calibrators').splitlines()
                          if field in fields]
    logger.info("Secondary calibrators: {0}".format(secondary_cals))
    #
    # Get flux calibrator fields if they are not in config
    #
    if config.get('Calibrators', 'Flux Calibrators') == '':
        flux_cals = []
        logger.info("Looking for flux calibrators...")
        with open('listobs.txt') as f:
            for line in f:
                for field in fields:
                    if field in line and "CALIBRATE_FLUX" in line \
                      and field not in flux_cals:
                        flux_cals.append(field)
        logger.info("Done")
    else:
        flux_cals = [field for field in config.get('Calibrators','Flux Calibrators').splitlines()
                     if field in fields]
    logger.info("Flux calibrators: {0}".format(flux_cals))
    #
    # Check that flux calibrators are in primary calibrator list
    # if not, add them
    #
    for flux_cal in flux_cals:
        if flux_cal not in primary_cals:
            primary_cals.append(flux_cal)
    #
    # Get science targets
    # 
    science_targets = []
    logger.info("Looking for science targets...")
    for field in fields:
        if field not in primary_cals+secondary_cals:
            science_targets.append(field)
    logger.info("Done")
    logger.info("Science targets: {0}".format(science_targets))
    #
    # create directories for figures
    #
    if not os.path.isdir('calibrator_figures'):
        logger.info("Creating calibrator_figures directory...")
        os.makedirs('calibrator_figures')
        logger.info("Done.")
    if not os.path.isdir('science_figures'):
        logger.info("Creating science_figures directory...")
        os.makedirs('science_figures')
        logger.info("Done.")
    if not os.path.isdir('plotcal_figures'):
        logger.info("Creating plotcal_figures directory...")
        os.makedirs('plotcal_figures')
        logger.info("Done.")
    return (my_cont_spws, my_line_spws, flux_cals, primary_cals,
            secondary_cals, science_targets, refant)

def preliminary_flagging(vis='', my_line_spws='', my_cont_spws='',
                         shadow_tolerance=0.0, quackinterval=10.0,
                         config=None):
    """
    Perform preliminary flagging: shadowed antennas, quack,
    flags from configuration file, interpolations from configuration
    file, then extend flags as necessary.

    Inputs: vis, my_line_spws, my_cont_spws, config
      vis :: string
        The measurement set
      my_line_spws :: string
        comma-separated string of line spws
      my_cont_spws :: string
        comma-separated string of continuum spws
      shadow_tolerance :: scalar
        The overlap tolerance used for shadow flagging. Flag
        any data with projected baseline separation less than
        r_1 + r_2 - shadow_tolerance
        where r_1 and r_2 are the radii of the antennas.
      quackinterval :: scalar
        The amount of time in seconds to flag at the beginning
        of each scan.
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None")
    #
    # save initial flag state
    #
    logger.info("Saving flag state...")
    casa.flagmanager(vis=vis, mode='save',
                     versionname='starting_flags_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    logger.info("Done.")
    #
    # Flag shadowed antennas
    #
    logger.info("Flagging shadowed antennas...")
    casa.flagdata(vis=vis, mode='shadow', tolerance=shadow_tolerance,
                  flagbackup=False, extendflags=False)
    logger.info("Done.")
    #
    # Flag the beginning of each scan
    #
    logger.info("Flagging the beginning of each scan (quack)...")
    casa.flagdata(vis=vis, mode='quack', quackinterval=quackinterval,
                  flagbackup=False, extendflags=False)
    logger.info("Done.")
    #
    # Flag scans from configuration file
    #
    scan = config.get("Flags", "Scan")
    if antenna != '':
        logger.info("Flagging scans from configuration file: {0}".format(scan))
        casa.flagdata(vis=vis, mode='manual', scan=scan,
                      flagbackup=False, extendflags=False)
        logger.info("Done.")
    #
    # Flag antennas from configuration file
    #
    antenna = config.get("Flags", "Antenna")
    if antenna != '':
        logger.info("Flagging antennas from configuration file: {0}".format(antenna))
        casa.flagdata(vis=vis, mode='manual', antenna=antenna,
                      flagbackup=False, extendflags=False)
        logger.info("Done.")
    #
    # Flag line channels from configuration file
    #
    badchans = config.get("Flags", "Line Channels").split(',')
    if badchans[0] != '':
        logger.info("Flagging line channels from configuration file: {0}".format(badchans))
        badchans = ';'.join(badchans)
        line_spws = ','.join([i+':'+badchans for i in my_line_spws.split(',')])
        casa.flagdata(vis=vis, mode='manual', spw=line_spws,
                      flagbackup=False, extendflags=False)
        logger.info("Done.")
    #
    # Flag continuum channels from configuration file
    #
    badchans = config.get("Flags", "Continuum Channels").split(',')
    if badchans[0] != '':
        logger.info("Flagging continuum channels from configuration file: {0}".format(badchans))
        badchans = ';'.join(badchans)
        cont_spws = ','.join([i+':'+badchans for i in my_cont_spws.split(',')])
        casa.flagdata(vis=vis, mode='manual', spw=cont_spws,
                      flagbackup=False, extendflags=False)
        logger.info("Done.")
    #
    # Interpolate through bad line channels
    #
    badchans = config.get("Interpolate", "Line Channels").split(',')
    if badchans[0] != '':
        logger.info("Interpolating through line channels from configuration file: {0}".format(badchans))
        line_spws = [int(i) for i in my_line_spws.split(',')]
        badchans = np.array([int(i) for i in badchans])
        for line_spw in line_spws:
            casa.ms.open(vis, nomodify=False)
            logger.info("Working on spw {0}".format(line_spw))
            casa.ms.selectinit(datadescid=line_spw)
            foo = casa.ms.getdata(['data'])
            #
            # Loop over time and polarization axes
            #
            for time_ind in range(len(foo['data'][0,0,:])):
                for pol_ind in range(len(foo['data'][:,0,0])):
                    # calculate amplitude from complex visibilities
                    amp = np.abs(foo['data'][pol_ind,:,time_ind])
                    # interpolate
                    amp[badchans] = np.interp(badchans, np.delete(range(len(amp)), badchans),
                                              np.delete(amp, badchans))
                    # calculate and unwrap phase from complex visibility
                    phase = np.unwrap(np.angle(foo['data'][pol_ind,:,time_ind]))
                    # interpolate
                    phase[badchans] = np.interp(badchans, np.delete(range(len(phase)), badchans),
                                                np.delete(phase, badchans))
                    # re-calculate complex visibilities
                    foo['data'][pol_ind,:,time_ind] = amp*np.cos(phase) + 1.j*amp*np.sin(phase)
            # save new data
            casa.ms.putdata(foo)
            casa.ms.done()
            foo = 0
            gc.collect()
    #
    # Interpolate through bad continuum channels
    #
    badchans = config.get("Interpolate", "Continuum Channels").split(',')
    if badchans[0] != '':
        logger.info("Interpolating through continuum channels from configuration file: {0}".format(badchans))
        cont_spws = [int(i) for i in my_cont_spws.split(',')]
        badchans = np.array([int(i) for i in badchans])
        for cont_spw in cont_spws:
            casa.ms.open(vis, nomodify=False)
            logger.info("Working on spw {0}".format(cont_spw))
            casa.ms.selectinit(datadescid=cont_spw)
            foo = casa.ms.getdata(['data'])
            #
            # Loop over time and polarization axes
            #
            for time_ind in range(len(foo['data'][0,0,:])):
                for pol_ind in range(len(foo['data'][:,0,0])):
                    # calculate amplitude from complex visibilities
                    amp = np.abs(foo['data'][pol_ind,:,time_ind])
                    # interpolate
                    amp[badchans] = np.interp(badchans, np.delete(range(len(amp)), badchans),
                                              np.delete(amp, badchans))
                    # calculate and unwrap phase from complex visibility
                    phase = np.unwrap(np.angle(foo['data'][pol_ind,:,time_ind]))
                    # interpolate
                    phase[badchans] = np.interp(badchans, np.delete(range(len(phase)), badchans),
                                                np.delete(phase, badchans))
                    # re-calculate complex visibilities
                    foo['data'][pol_ind,:,time_ind] = amp*np.cos(phase) + 1.j*amp*np.sin(phase)
            # save new data
            casa.ms.putdata(foo)
            casa.ms.done()
            foo = 0
            gc.collect()
    #
    # Extend the flags
    #
    logger.info("Extending flags...")
    casa.flagdata(vis=vis, mode='extend', extendpols=True,
                  growtime=80.0, growfreq=50.0, growaround=True,
                  flagbackup=False)
    logger.info("Done.")
    
    #
    # Save the flags
    #
    logger.info("Saving flag state...")
    casa.flagmanager(vis=vis, mode='save',
                     versionname='preliminary_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    logger.info("Done.")

def auto_flag_calibrators(vis='', primary_cals=[], secondary_cals=[]):
    """
    Perform automatic flagging of calibrators using rflag on
    calibrated data, then extend the flags as necessary.

    Inputs: vis, primary_cals, secondary_cals
      vis :: string
        The measurement set
      primary_cals :: list of strings
        List of primary calibrator names (including flux calibrators)
      secondary_cals :: list of strings
        List of secondary calibrator names

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check if calibrators have corrected datacolumn
    #
    field = ','.join(primary_cals+secondary_cals)
    stat = None
    logger.info("Checking if ms contains corrected data column...")
    stat = casa.visstat(vis=vis, field=field, spw='0', datacolumn='corrected')
    if stat is None:
        logger.critical("Done. ms does not contain corrected data column. Skipping.")
        return
    logger.info("Done. ms does contain corrected data column.")
    #
    # Run rflag on calibrated data
    #
    logger.info("Running rflag on corrected data column...")
    casa.flagdata(vis=vis, mode='rflag', field=field,
                  flagbackup=False, datacolumn='corrected',
                  extendflags=False)
    logger.info("Done.")
    #
    # Extend the flags
    #
    logger.info("Extending flags...")
    casa.flagdata(vis=vis, mode='extend', extendpols=True,
                  growtime=80.0, growfreq=50.0, growaround=True,
                  flagbackup=False)
    logger.info("Done.")
    #
    # Save the flags
    #
    logger.info("Saving flag state...")
    casa.flagmanager(vis=vis, mode='save',
                     versionname='autoflag_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    logger.info("Done.")

def gen_calibrator_plots(vis='',primary_cals=[],secondary_cals=[],
                         config=None):
    """
    Generate diagnostic visibility plots for calibrators

    Inputs: vis, primary_cals, secondary_cals, config
      vis :: string
        The measurement set
      primary_cals :: list of strings
        list of primary calibrator names (must include flux calibrators)
      secondary_cals :: list of strings
        list of secondary calibrator names
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None")
    #
    # check if calibrators have corrected datacolumn
    #
    field = ','.join(primary_cals+secondary_cals)
    stat = None
    logger.info("Checking if ms contains corrected data column...")
    stat = casa.visstat(vis=vis, spw='0', field=field, datacolumn='corrected')
    if stat is None:
        logger.info("Done. ms does not contain corrected data column.")
        datacolumn='data'
    else:
        logger.info("Done. ms does contain corrected data column.")
        datacolumn='corrected'
    #
    # Generate the plots
    #
    logger.info("Generating plots for manual inspection...")
    plotnum=0
    plots = []
    for field in primary_cals+secondary_cals:
        #
        # Phase vs. Amplitude
        #
        casa.plotms(vis=vis, xaxis='amp', yaxis='phase', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    plotfile='calibrator_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'amp', 'yaxis':'phase', 'avgtime':'', 'avgchannel':''})
        plotnum += 1
        #
        # Amplitude vs UV-distance (in wavelength units)
        #
        casa.plotms(vis=vis, xaxis='uvwave', yaxis='amp', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    plotfile='calibrator_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'uvwave', 'yaxis':'amp', 'avgtime':'', 'avgchannel':''})
        plotnum += 1
         #
        # Amplitude vs Time
        #
        casa.plotms(vis=vis, xaxis='time', yaxis='amp', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    avgchannel='1e7', 
                    plotfile='calibrator_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'time', 'yaxis':'amp', 'avgtime':'', 'avgchannel':'1e7'})
        plotnum += 1
        #
        # Phase vs Time
        #
        casa.plotms(vis=vis, xaxis='time', yaxis='phase', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    avgchannel='1e7', 
                    plotfile='calibrator_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'time', 'yaxis':'phase', 'avgtime':'', 'avgchannel':'1e7'})
        plotnum += 1
        #
        # Amplitude vs Channel
        #
        casa.plotms(vis=vis, xaxis='channel', yaxis='amp', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    avgtime='1e7', 
                    plotfile='calibrator_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'channel', 'yaxis':'amp', 'avgtime':'1e7', 'avgchannel':''})
        plotnum += 1
        #
        # Phase vs Channel
        # 
        casa.plotms(vis=vis, xaxis='channel', yaxis='phase', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    avgtime='1e7', 
                    plotfile='calibrator_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'channel', 'yaxis':'phase', 'avgtime':'1e7', 'avgchannel':''})
        plotnum += 1
    logger.info("Done.")
    #
    # Generate PDF to display plots
    #
    logger.info("Generating PDF...")
    num_plots = plotnum
    iplot = 0
    with open('calibrator_figures.tex', 'w') as f:
        f.write(r"\documentclass{article}"+"\n")
        f.write(r"\usepackage{graphicx}"+"\n")
        f.write(r"\usepackage[margin=0.1cm]{geometry}"+"\n")
        f.write(r"\begin{document}"+"\n")
        f.write(r"\begin{figure}"+"\n")
        f.write(r"\centering"+"\n")
        for plotnum in range(num_plots):
            fnames=glob.glob("calibrator_figures/{0}_*.png".format(plotnum))
            fnames = natural_sort(fnames)
            for fname in fnames:
                if iplot > 0 and iplot % 6 == 0:
                    f.write(r"\end{figure}"+"\n")
                    f.write(r"\clearpage"+"\n")
                    f.write(r"\begin{figure}"+"\n")
                    f.write(r"\centering"+"\n")
                elif iplot > 0 and iplot % 2 == 0:
                    f.write(r"\end{figure}"+"\n")
                    f.write(r"\begin{figure}"+"\n")
                    f.write(r"\centering"+"\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+fname+"}\n")
                iplot+=1
        f.write(r"\end{figure}"+"\n")
        f.write(r"\end{document}"+"\n")
    os.system('pdflatex -interaction=batchmode calibrator_figures.tex')
    logger.info("Done.")
    #
    # Save plot list to a pickle object
    #
    logger.info("Saving plot list to pickle...")
    with open('calibrator_plots.pkl', 'w') as f:
        pickle.dump(plots, f)
    logger.info("Done.")

def flag(vis='', all_fields=[]):
    """
    Interactively flag some data

    Inputs: vis, all_fields
      vis :: string
        The measurement set
      all_fields :: list of strings
        List of field names to be flagged

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Build list of flag commands
    #
    flag_commands = []
    while True:
        #
        # Prompt user for field, scan, spectral window, time range,
        # antenna, and correlation to flag
        #
        print("Field? Empty = {0}".format(','.join(all_fields)))
        field = raw_input()
        if field == '':
            field = ','.join(all_fields)
        if field not in ','.join(all_fields):
            print('{0} is not in {1}'.format(field, all_fields))
            continue
        print("Scan? Empty = all scans (ex. 0 to flag scan 0, 1~3 to flag scans 1, 2, and 3)")
        scan = raw_input()
        print("Spectral window and channels? Empty = all spws (ex. 2:100~150;200:250,5:1200~1210 is spw 2, chans 100 to 150 and 200 to 250, and spw 5 chans 1200 to 1210)")
        spw = raw_input()
        print("Time range? Empty = all times (ex. 10:23:45~10:23:55)")
        timerange = raw_input()
        print("Antenna or baseline? Empty = all antennas/baselines (ex. CA01 to flag ant 1 or CA01&CA02 to flag 1-2 baseline)")
        antenna = raw_input()
        print("Correlation? Empty = all correlations (XX,YY,XY,YX)")
        correlation = raw_input()
        #
        # Build flag command
        #
        parts = []
        if len(field) > 0:
            parts.append("field='{0}'".format(field))
        if len(scan) > 0:
            parts.append("scan='{0}'".format(scan))
        if len(spw) > 0:
            parts.append("spw='{0}'".format(spw))
        if len(timerange) > 0:
            parts.append("timerange='{0}'".format(timerange))
        if len(antenna) > 0:
            parts.append("antenna='{0}'".format(antenna))
        if len(correlation) > 0:
            parts.append("correlation='{0}'".format(correlation))
        flag_commands.append(' '.join(parts))
        #
        # Confirm with user, or append more flag commands
        #
        print("Will execute:")
        print("flagdata(vis='{0}',mode='list',flagbackup=False,extendflags=False,".format(vis))
        for icmd,cmd in enumerate(flag_commands):
            if len(flag_commands) == 1:
                print("         inpfile=[\"{0}\"])".format(cmd))
            elif icmd == 0:
                print("         inpfile=[\"{0}\",".format(cmd))
            elif icmd == len(flag_commands)-1:
                print("                  \"{0}\"])".format(cmd))
            else:
                print("                  \"{0}\",".format(cmd))
        print("Proceed [y/n] or add another flag command [a]?")
        go = raw_input()
        #
        # Execute flag command
        #
        if go.lower() == 'y':
            logger.info("Executing:")
            logger.info("flagdata(vis='{0}',mode='list',flagbackup=False,extendflags=False,".format(vis))
            for icmd,cmd in enumerate(flag_commands):
                if len(flag_commands) == 1:
                    logger.info("         inpfile=[\"{0}\"])".format(cmd))
                elif icmd == 0:
                    logger.info("         inpfile=[\"{0}\",".format(cmd))
                elif icmd == len(flag_commands)-1:
                    logger.info("                  \"{0}\"])".format(cmd))
                else:
                    logger.info("                  \"{0}\",".format(cmd))
            #
            # Save flag command to manual flags list
            #
            with open('manual_flags.txt','a') as f:
                for cmd in flag_commands:
                    f.write(time.strftime('%Y%m%d%H%M%S', time.gmtime())+': '+cmd+'\n')
            casa.flagdata(vis=vis, mode='list', flagbackup=False, extendflags=False,
                          inpfile=flag_commands)
            break
        #
        # Append another flag command
        #
        elif go.lower() == 'a':
            continue
        #
        # Quit flagging
        #
        else:
            print("Aborting...")
            break

def manual_flag_calibrators(vis='', primary_cals=[], secondary_cals=[],
                            science_targets=[], config=None):
    """
    Interactively plot and flag the calibrators

    Inputs: vis, primary_cals, secondary_cals, science_targets, config
      vis :: string
        The measurement set
      primary_cals :: list of strings
        list of primary calibrator names (must include flux calibrators)
      secondary_cals :: list of strings
        list of secondary calibrator names
      science_targets :: list of strings
        list of science target names
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None")
    #
    # check if calibrators have corrected datacolumn
    #
    field = ','.join(primary_cals+secondary_cals)
    stat = None
    logger.info("Checking if ms contains corrected data column...")
    stat = casa.visstat(vis=vis, spw='0', field=field, datacolumn='corrected')
    if stat is None:
        logger.info("Done. ms does not contain corrected data column.")
        datacolumn='data'
    else:
        logger.info("Done. ms does contain corrected data column.")
        datacolumn='corrected'
    #
    # Read the plot list from the pickle object
    #
    logger.info("Reading plot list from pickle...")
    with open('calibrator_plots.pkl','r') as f:
        plots = pickle.load(f)
    num_plots = len(plots)
    logger.info("Done.")
    #
    # Display menu option to user
    #
    logger.info("Please inspect calibrator_plots.pdf then perform manual flagging.")
    while True:
        print("f - flag some data")
        print("plot id number - generate interactive version of plot with this id")
        print("quit - end this flagging session")
        answer = raw_input()
        #
        # Flag some data
        #
        if answer.lower() == 'f':
            # if we are going to flag something in all fields for the
            # calibrators, we might as well flag it in the science targets
            # too
            flag(vis,all_fields=primary_cals+secondary_cals+science_targets)
        #
        # Stop interactively plotting and flagging
        #
        elif answer.lower() == 'quit':
            break
        #
        # Generate interactive plot
        #
        else:
            try:
                plotid = int(answer)
            except ValueError:
                print("Plot ID not valid!")
                continue
            if plotid >= num_plots:
                print("Plot ID not valid!")
                continue
            casa.plotms(vis=vis, xaxis=plots[plotid]['xaxis'], yaxis=plots[plotid]['yaxis'],
                        field=plots[plotid]['field'], ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=config.get('Polarization','Polarization'),
                        title='PlotID: {0} Field: {1}'.format(plotid,plots[plotid]['field']),
                        avgchannel=plots[plotid]['avgchannel'],
                        avgtime=plots[plotid]['avgtime'])
    #
    # Save the flags
    #
    logger.info("Saving flag state...")
    casa.flagmanager(vis=vis, mode='save',
                     versionname='manualflag_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    logger.info("Done.")

def calibrate_calibrators(vis='', primary_cals=[], secondary_cals=[],
                          flux_cals=[], my_line_spws='',
                          my_cont_spws='', refant='', calwt=True,
                          antpos=True, gaincurve=True, opacity=True,
                          config=None):
    """
    Calculate calibration solutions (bandpass, delays, complex gains)
    and apply the calibration solutions to the calibrators

    Inputs: vis, primary_cals, secondary_cals, flux_cals,
            my_line_spws, my_cont_spws, refant, calwt, config
      vis :: string
        The measurement set
      primary_cals :: list of strings
        list of primary calibrator names (must include flux calibrators)
      secondary_cals :: list of strings
        list of secondary calibrator names
      flux_cals :: list of strings
        list of flux calibrator names
      my_line_spws :: string
        comma-separated string of line spectral windows
      my_cont_spws :: string
        comma-separated string of continuum spectral windows
      refant :: string
        reference antenna
      calwt :: boolean
        if True, apply calibration weights to data
      antpos :: boolean
        if True, compute antenna position corrections (only for VLA)
      gaincurve :: boolean
        if True, compute gain curve and antenna efficiency corrections
      opacity :: boolean
        if True, compute opacity corrections
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: caltables, calfields
      caltables :: list of strings
        The names of the antpos, gaincurve, and/or opacity calibration
        tables if computed
      calfields :: list of strings
        The field numbers for the antpos, gaincurve, and/or opacity
        calibration tables if computed. These are always empty
        strings because they aren't associated with a single field.
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None")
    #
    # set the model for the flux calibrators
    #
    logger.info("Setting the flux calibrator models...")
    for flux_cal in flux_cals:
        #
        # If flux calibrator model is supplied in config, use that
        #
        manual_flux_cals = config.get('Flux Calibrator Models','Name').splitlines()
        if flux_cal in manual_flux_cals:
            # get index of flux_cal in manual_flux_cals
            flux_idx = manual_flux_cals.index(flux_cal)
            # get reference frequency
            reffreq = config.get('Flux Calibrator Models','Reference Frequency').splitlines()
            reffreq = reffreq[flux_idx]
            # get fluxdensity and convert to proper units
            fluxdensity = config.get('Flux Calibrator Models','Log Flux Density').splitlines()
            fluxdensity = [10.**float(fluxdensity[flux_idx]),0.,0.,0.]
            # get spectral index coefficients
            spix = config.get('Flux Calibrator Models','Spectral Index Coefficients').splitlines()
            spix = spix[flux_idx]
            spix = [float(i) for i in spix.split(',')]
            # Run setjy in manual mode
            casa.setjy(vis=vis, field=flux_cal, scalebychan=True,
                       standard='manual', fluxdensity=fluxdensity,
                       spix=spix, reffreq=reffreq)
        #
        # Otherwise, use CASA model
        #
        else:
            casa.setjy(vis=vis, field=flux_cal, scalebychan=True)
    logger.info("Done.")
    caltables = []
    calfields = []
    #
    # Correct antenna positions
    #
    if antpos:
        logger.info("Calculating antenna position correction...")
        if os.path.isdir('antpos.cal'):
            casa.rmtables('antpos.cal')
        casa.gencal(vis=vis, caltable='antpos.cal', caltype='antpos')
        if os.path.isdir('antpos.cal'):
            caltables += ['antpos.cal']
            calfields += ['']
        logger.info("Done.")
    #
    # Correct for gaincurve and antenna efficiencies
    #
    if gaincurve:
        logger.info("Calculating gain curve and antenna efficiencies...")
        if os.path.isdir('gaincurve.cal'):
            casa.rmtables('gaincurve.cal')
        casa.gencal(vis=vis, caltable='gaincurve.cal', caltype='gceff')
        caltables += ['gaincurve.cal']
        calfields += ['']
        logger.info("Done.")
    #
    # Correct for atmospheric opacity
    #
    if opacity:
        logger.info("Calculating opacity correction...")
        if os.path.isdir('opacity.cal'):
            casa.rmtables('opacity.cal')
        myTau = casa.plotweather(vis=vis,doPlot=False)
        allspws = ','.join([str(spw) for spw in range(len(myTau))])
        casa.gencal(vis=vis, caltable='opacity.cal', caltype='opac',
                    parameter=myTau, spw=allspws)
        caltables += ['opacity.cal']
        calfields += ['']
        logger.info("Done.")
    #
    # pre-bandpass calibration delay calibration on primary calibrators
    # (linear slope in phase vs frequency)
    #
    field = ','.join(primary_cals)
    logger.info("Calculating delay calibration table for primary calibrators...")
    if os.path.isdir('delays.Kcal'):
        casa.rmtables('delays.Kcal')
    casa.gaincal(vis=vis, caltable='delays.Kcal', field=field,
                 refant=refant, gaintype='K', minblperant=1,
                 gaintable = caltables)
    if not os.path.isdir('delays.Kcal'):
        logger.critical('Problem with delay calibration')
        raise ValueError('Problem with delay calibration!')
    logger.info("Done.")
    #
    # integration timescale phase calibration
    # (phase vs time)
    #
    logger.info("Calculating phase calibration table on integration timescales for primary calibrators...")
    if os.path.isdir('phase_int.Gcal0'):
        casa.rmtables('phase_int.Gcal0')
    casa.gaincal(vis=vis, caltable="phase_int.Gcal0", field=field,
                 solint="int", calmode="p", refant=refant,
                 gaintype="G", minsnr=2.0, minblperant=1,
                 gaintable=caltables + ['delays.Kcal'])
    if not os.path.isdir('phase_int.Gcal0'):
        logger.critical('Problem with integration-timescale phase calibration')
        raise ValueError('Problem with integration-timescale phase calibration!')
    logger.info("Done.")
    #
    # Generate phase_int.Gcal0 plots
    #
    casa.plotms(vis='phase_int.Gcal0', xaxis='time', yaxis='phase',
                field=field, iteraxis='spw',
                coloraxis='antenna',
                title='Phase GCal (integration; pre-bandpass)',
                plotfile='plotcal_figures/0_phase_int.Gcal0.png',
                overwrite=True, showgui=False, exprange='all')
    #
    # bandpass calibration for continuum spws. Combine all scans,
    # average some channels as defined in configuration file
    #
    logger.info("Calculating bandpass calibration table for primary calibrators...")
    if os.path.isdir('bandpass.Bcal'):
        casa.rmtables('bandpass.Bcal')
    chan_avg = config.get('Bandpass Channel Average','Continuum Channels')
    if chan_avg == '':
        solint='inf'
    else:
        solint='inf,{0}chan'.format(chan_avg)
    casa.bandpass(vis=vis, caltable='bandpass.Bcal', field=field,
                  spw=my_cont_spws, refant=refant, solint=solint,
                  combine='scan', solnorm=True, minblperant=1,
                  gaintable=caltables + ['delays.Kcal', 'phase_int.Gcal0'])
    #
    # bandpass calibration for line spws. Combine all scans,
    # average some channels as defined in configuration file,
    # append to continuum channel bandpass calibration table
    #
    chan_avg = config.get('Bandpass Channel Average', 'Line Channels')
    if chan_avg == '':
        solint='inf'
    else:
        solint='inf,{0}chan'.format(chan_avg)
    casa.bandpass(vis=vis, caltable='bandpass.Bcal', field=field,
                  spw=my_line_spws, refant=refant, solint=solint,
                  combine='scan', solnorm=True, minblperant=1, append=True,
                  gaintable=caltables + ['delays.Kcal', 'phase_int.Gcal0'])
    if not os.path.isdir('bandpass.Bcal'):
        logger.critical('Problem with bandpass calibration')
        raise ValueError('Problem with bandpass calibration!')
    logger.info("Done.")
    #
    # Generate bandpass.Bcal plots
    #
    casa.plotms(vis='bandpass.Bcal', xaxis='channel', yaxis='amplitude',
                field=field, iteraxis='spw',
                coloraxis='antenna',
                title='Bandpass BCal',
                plotfile='plotcal_figures/1_bandpass.Bcal.png',
                overwrite=True, showgui=False, exprange='all')
    #
    # integration timescale phase corrections for all calibrators
    # required for accurate amplitude calibration
    #
    field = ','.join(primary_cals+secondary_cals)
    logger.info("Re-calculating the phase calibration table on integration timescales for all calibrators...")
    if os.path.isdir('phase_int.Gcal1'):
        casa.rmtables('phase_int.Gcal1')
    casa.gaincal(vis=vis, caltable="phase_int.Gcal1", field=field,
                solint="int", calmode="p", refant=refant,
                gaintype="G", minsnr=2.0, minblperant=1,
                gaintable=caltables + ['delays.Kcal', 'bandpass.Bcal'])
    if not os.path.isdir('phase_int.Gcal1'):
        logger.critical('Problem with integration-timescale phase calibration')
        raise ValueError('Problem with integration-timescale phase calibration!')
    logger.info("Done.")
    #
    # Generate phase_int.Gcal1 plots
    #
    casa.plotms(vis='phase_int.Gcal1', xaxis='time', yaxis='phase',
                field=field, iteraxis='spw',
                coloraxis='antenna',
                title='Phase GCal (integration; post-bandpass)',
                plotfile='plotcal_figures/2_phase_int.Gcal1.png',
                overwrite=True, showgui=False, exprange='all')
    #
    # scan timescale phase corrections for all calibrators
    # required to apply to science targets
    #
    logger.info("Calculating the phase calibration table on scan timescales for all calibrators...")
    if os.path.isdir('phase_scan.Gcal'):
        casa.rmtables('phase_scan.Gcal')
    casa.gaincal(vis=vis, caltable="phase_scan.Gcal", field=field,
                 solint="inf", calmode="p", refant=refant,
                 gaintype="G", minsnr=2.0, minblperant=1,
                 gaintable=caltables + ['delays.Kcal', 'bandpass.Bcal'])
    if not os.path.isdir('phase_scan.Gcal'):
        logger.critical('Problem with scan-timescale phase calibration')
        raise ValueError('Problem with scan-timescale phase calibration!')
    logger.info("Done.")
    #
    # Generate phase_scan.Gcal plots
    #
    casa.plotms(vis='phase_scan.Gcal', xaxis='time', yaxis='phase',
                field=field, iteraxis='spw',
                coloraxis='antenna',
                title='Phase GCal (scan; post-bandpass)',
                plotfile='plotcal_figures/3_phase_scan.Gcal.png',
                overwrite=True, showgui=False, exprange='all')
    #
    # scan timescale amplitude corrections using
    # integration timescale phase calibration
    #
    logger.info("Calculating the amplitude calibration table on scan timescales for all calibrators...")
    if os.path.isdir('apcal_scan.Gcal'):
        casa.rmtables('apcal_scan.Gcal')
    casa.gaincal(vis=vis, caltable="apcal_scan.Gcal", field=field,
                 solint="inf", calmode="ap", refant=refant,
                 minsnr=2.0, minblperant=1,
                 gaintable=caltables + ['delays.Kcal', 'bandpass.Bcal', 'phase_int.Gcal1'])
    if not os.path.isdir('apcal_scan.Gcal'):
        logger.critical('Problem with amplitude calibration')
        raise ValueError('Problem with amplitude calibration!')
    logger.info("Done.")
    #
    # Generate apcal_scan.Gcal plots
    #
    casa.plotms(vis='apcal_scan.Gcal', xaxis='time', yaxis='amplitude',
                field=field, iteraxis='spw',
                coloraxis='antenna',
                title='Amplitude GCal (scan; post-bandpass)',
                plotfile='plotcal_figures/4_apcal_scan.Gcal.png',
                overwrite=True, showgui=False, exprange='all')
    #
    # set the flux scale
    #
    logger.info("Calculating the flux calibration table...")
    if os.path.isdir('flux.cal'):
        casa.rmtables('flux.cal')
    casa.fluxscale(vis=vis, caltable="apcal_scan.Gcal", fluxtable="flux.cal",
                   reference=','.join(flux_cals), incremental=True)
    if not os.path.isdir('flux.cal'):
        logger.critical('Problem with flux calibration')
        raise ValueError('Problem with flux calibration!')
    logger.info("Done.")
    #
    # Generate flux.cal plots
    #
    casa.plotms(vis='flux.cal', xaxis='channel', yaxis='amplitude',
                field=field, iteraxis='spw',
                coloraxis='antenna',
                title='Flux Cal',
                plotfile='plotcal_figures/5_flux.cal.png',
                overwrite=True, showgui=False, exprange='all')
    #
    # apply calibration solutions to calibrators
    #
    logger.info("Applying calibration tables to all calibrators...")
    for field in primary_cals+secondary_cals:
        casa.applycal(vis=vis, field=field, calwt=calwt,
                      gaintable=caltables + ['delays.Kcal', 'bandpass.Bcal',
                                             'phase_int.Gcal1', 'apcal_scan.Gcal',
                                             'flux.cal'],
                      gainfield=calfields + ['', '', field, field, field],
                      flagbackup=False)
    logger.info("Done.")
    #
    # save flag state
    #
    casa.flagmanager(vis=vis, mode='save',
                     versionname='calibrate_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    #
    # Generate PDF of plotcal figures
    #
    logger.info("Generating tex document...")
    iplot = 0
    with open('plotcal_figures.tex','w') as f:
        f.write(r"\documentclass{article}"+"\n")
        f.write(r"\usepackage{graphicx,subfig}"+"\n")
        f.write(r"\usepackage[margin=0.1cm]{geometry}"+"\n")
        f.write(r"\begin{document}"+"\n")
        f.write(r"\begin{figure}"+"\n")
        f.write(r"\centering"+"\n")
        fnames = glob.glob("plotcal_figures/{0}_*.png".format(plotnum))
        fnames = natural_sort(fnames)
        for fname in fnames:
            if iplot > 0 and iplot % 6 == 0:
                f.write(r"\end{figure}"+"\n")
                f.write(r"\clearpage"+"\n")
                f.write(r"\begin{figure}"+"\n")
                f.write(r"\centering"+"\n")
            elif iplot > 0 and iplot % 2 == 0:
                f.write(r"\end{figure}"+"\n")
                f.write(r"\begin{figure}"+"\n")
                f.write(r"\centering"+"\n")
            f.write(r"\subfloat[\texttt{"+fname.replace("_",r"\_")+"}]{")
            f.write(r"\includegraphics[width=0.45\textwidth]{"+fname+"}}"+"\n")
            iplot+=1
        f.write(r"\end{figure}"+"\n")
        f.write(r"\end{document}"+"\n")
    os.system('pdflatex -interaction=batchmode plotcal_figures.tex')
    logger.info("Done.")
    return caltables, calfields

def calibrate_sciencetargets(caltables, calfields, vis='',
                             science_targets=[], calwt=True):
    """
    Apply calibration solutions to science targets

    Inputs: caltables, calfields, vis, science_targets, calwt
      caltables :: list of strings
        Returned by calibrate_calibrators()
      calfields :: list of strings
        Returned by calibrate_calibrators()
      vis :: string
        The measurement set
      science_targets :: list of strings
        list of science target names
      calwt :: boolean
        if True, apply calibration weights to data
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    for field in science_targets:
        #
        # use all fields in delays and bandpass
        # use nearest field in complex gains and flux
        #
        logger.info("Applying calibration solutions to {0}".format(field))
        casa.applycal(vis=vis,field=field, calwt=calwt,
                      gaintable=caltables+['delays.Kcal', 'bandpass.Bcal',
                                           'phase_scan.Gcal', 'apcal_scan.Gcal',
                                           'flux.cal'],
                      gainfield=calfields+['','','nearest','nearest','nearest'],
                      flagbackup=False)
    #
    # save the flags
    #
    logger.info("Saving flag state...")
    casa.flagmanager(vis=vis, mode='save',
                     versionname='calibrate_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    logger.info("Done.")

def auto_flag_sciencetargets(vis='', science_targets=[]):
    """
    Perform automatic flagging of calibrated science targets using
    rflag

    Inputs: vis, science_targets
      vis :: string
        The measurement set
      science_targets :: list of strings
        List of science target names

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Perform automatic flagging on science targets
    #
    field = ','.join(science_targets)
    datacolumn='corrected'
    logger.info("Running rflag on all correlations...")
    casa.flagdata(vis=vis, mode='rflag', field=field,
                  flagbackup=False, datacolumn=datacolumn,
                  extendflags=False)
    logger.info("Done.")
    #
    # Save the flags
    #
    logger.info("Saving flag state...")
    casa.flagmanager(vis=vis, mode='save',
                     versionname='autoflag_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    logger.info("Done.")

def gen_sciencetarget_plots(vis='', science_targets=[], config=None):
    """
    Generate science target visibility plots

    Inputs: vis, science_targets, config
      vis :: string
        The measurement set
      science_targets :: list of strings
        list of science target names
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    datacolumn='corrected'
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None")
    #
    # Generate the plots
    #
    logger.info("Generating plots for manual inspection...")
    plotnum=0
    plots = []
    for field in science_targets:
        #
        # Amplitude vs UV-distance (in wavelength units)
        #
        casa.plotms(vis=vis, xaxis='uvwave', yaxis='amp', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    plotfile='science_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'uvwave', 'yaxis':'amp', 'avgtime':'', 'avgchannel':''})
        plotnum += 1
        #
        # Amplitude vs Time
        #
        casa.plotms(vis=vis, xaxis='time', yaxis='amp', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    avgchannel='1e7', 
                    plotfile='science_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'time', 'yaxis':'amp', 'avgtime':'', 'avgchannel':'1e7'})
        plotnum += 1
        #
        # Amplitude vs Channel
        #
        casa.plotms(vis=vis, xaxis='channel', yaxis='amp', field=field, 
                    ydatacolumn=datacolumn, iteraxis='spw', 
                    coloraxis='baseline', correlation=config.get('Polarization', 'Polarization'), 
                    title='PlotID: {0} Field: {1}'.format(plotnum, field), 
                    avgtime='1e7', 
                    plotfile='science_figures/{0}.png'.format(plotnum), 
                    overwrite=True, showgui=False, exprange='all')
        plots.append({'field':field, 'xaxis':'channel', 'yaxis':'amp', 'avgtime':'1e7', 'avgchannel':''})
        plotnum += 1
    logger.info("Done.")
    #
    # Generate PDF to display plots
    #
    logger.info("Generating tex document...")
    num_plots = plotnum
    iplot = 0
    with open('science_figures.tex', 'w') as f:
        f.write(r"\documentclass{article}"+"\n")
        f.write(r"\usepackage{graphicx}"+"\n")
        f.write(r"\usepackage[margin=0.1cm]{geometry}"+"\n")
        f.write(r"\begin{document}"+"\n")
        f.write(r"\begin{figure}"+"\n")
        f.write(r"\centering"+"\n")
        for plotnum in range(num_plots):
            fnames=glob.glob("science_figures/{0}_*.png".format(plotnum))
            fnames = natural_sort(fnames)
            for fname in fnames:
                if iplot > 0 and iplot % 6 == 0:
                    f.write(r"\end{figure}"+"\n")
                    f.write(r"\clearpage"+"\n")
                    f.write(r"\begin{figure}"+"\n")
                    f.write(r"\centering"+"\n")
                elif iplot > 0 and iplot % 2 == 0:
                    f.write(r"\end{figure}"+"\n")
                    f.write(r"\begin{figure}"+"\n")
                    f.write(r"\centering"+"\n")
                f.write(r"\includegraphics[width=0.45\textwidth]{"+fname+"}\n")
                iplot+=1
        f.write(r"\end{figure}"+"\n")
        f.write(r"\end{document}"+"\n")
    os.system('pdflatex -interaction=batchmode science_figures.tex')
    logger.info("Done.")
    #
    # Save the plot list to the pickle object
    #
    logger.info("Saving plot list to pickle...")
    with open('science_plots.pkl', 'w') as f:
        pickle.dump(plots, f)
    logger.info("Done.")

def manual_flag_sciencetargets(vis='', science_targets=[], config=None):
    """
    Interactively plot and flag the science targets

    Inputs: vis, science_targets, config
      vis :: string
        The measurement set
      science_targets :: list of strings
        list of science target names
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    datacolumn='corrected'
    #
    # check config
    #
    if config is None:
        logger.critical("Error: Need to supply a config")
        raise ValueError("Config is None")
    #
    # Read plot list from pickle object
    #
    logger.info("Reading plot list from pickle...")
    with open("science_plots.pkl","r") as f:
        plots = pickle.load(f)
    num_plots = len(plots)
    logger.info("Done.")
    #
    # Prompt user with menu
    #
    logger.info("Please inspect science_plots.pdf then perform manual flagging.")
    while True:
        print("f - flag some data")
        print("plot id number - generate interactive version of plot with this id")
        print("quit - end this flagging session")
        answer = raw_input()
        #
        # Flag some data
        #
        if answer.lower() == 'f':
            flag(vis=vis, all_fields=science_targets)
        #
        # Stop flagging
        #
        elif answer.lower() == 'quit':
            break
        #
        # Generate plotms figure
        #
        else:
            try:
                plotid = int(answer)
            except ValueError:
                print("Plot ID not valid!")
                continue
            if plotid >= num_plots:
                print("Plot ID not valid!")
                continue
            casa.plotms(vis=vis, xaxis=plots[plotid]['xaxis'], yaxis=plots[plotid]['yaxis'],
                        field=plots[plotid]['field'], ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=config.get('Polarization','Polarization'),
                        title='PlotID: {0} Field: {1}'.format(plotid,plots[plotid]['field']),
                        avgchannel=plots[plotid]['avgchannel'],
                        avgtime=plots[plotid]['avgtime'])
    #
    # Save the flags
    #
    logger.info("Saving flag state...")
    casa.flagmanager(vis=vis, mode='save',
                     versionname='manualflag_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
    logger.info("Done.")

def split_fields(vis='', primary_cals=[], secondary_cals=[], science_targets=[]):
    """
    Split calibrated fields into own measurement sets with naming format:
    {field_name}_calibrated.ms

    Inputs: vis, primary_cals, secondary_cals, science_targets
      vis :: string
        The measurement set
      primary_cals :: list of strings
        list of primary calibrator names
      secondary_cals :: list of strings
        list of secondary calibrator names
      science_targets :: list of strings
        list of science target names

    Returns:
      Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    logger.info("Splitting fields...")
    for field in primary_cals+secondary_cals+science_targets:
        outputvis = '{0}_calibrated.ms'.format(field)
        logger.info("Splitting {0} to {1}".format(field,outputvis))
        # split and drop the flagged data
        casa.split(vis=vis, outputvis=outputvis, field=field, keepflags=False)
    logger.info("Done!")

def main(vis='', config_file='', shadow_tolerance=0.0, quackinterval=10.0,
         antpos=True, gaincurve=True, opacity=True, calwt=True, auto=''):
    """
    Run the calibration pipeline

    Inputs: vis, config_file, calwt, auto
      vis :: string
        The masurement set
      config_file :: string
        The filename of the configuration file for this project
      shadow_tolerance :: scalar
        The overlap tolerance used for shadow flagging. Flag
        any data with projected baseline separation less than
        r_1 + r_2 - shadow_tolerance
        where r_1 and r_2 are the radii of the antennas.
      quackinterval :: scalar
        The amount of time in seconds to flag at the beginning
        of each scan.
      antpos :: boolean
        if True, compute antenna position corrections (only for VLA)
      gaincurve :: boolean
        if True, compute gain curve and antenna efficiency corrections
      opacity :: boolean
        if True, compute opacity corrections
      calwt :: boolean
        if True, apply calibration weights to data
      auto :: string
        comma separated string of menu options to automatically
        perform

    Returns:
      Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Check inputs
    #
    if not os.path.isdir(vis):
        logger.critical('Measurement set not found!')
        raise ValueError('Measurement set not found!')
    if not os.path.exists(config_file):
        logger.critical('Configuration file not found')
        raise ValueError('Configuration file not found!')
    #
    # load configuration file
    #
    config = ConfigParser.ConfigParser()
    logger.info("Reading configuration file {0}".format(config_file))
    config.read(config_file)
    logger.info("Done.")
    #
    # initial setup
    #
    my_cont_spws,my_line_spws,flux_cals,primary_cals,secondary_cals,\
        science_targets,refant = setup(vis=vis,config=config)
    #
    # Prompt the user with a menu for each option, or auto-do them
    #
    auto_items = auto.split(',')
    auto_ind = 0
    caltables = []
    calfields = []
    while True:
        if len(auto) == 0:
            print("0. Preliminary flags (config file, etc.)")
            print("1. Auto-flag calibrator fields")
            print("2. Generate plotms figures for calibrator fields")
            print("3. Manually flag calibrator fields")
            print("4. Calculate and apply calibration solutions to calibrator fields")
            print("5. Apply calibration solutions to science fields")
            print("6. Auto-flag science fields")
            print("7. Generate plotms figures for science fields")
            print("8. Manually flag science fields")
            print("9. Split calibrated fields")
            print("q [quit]")
            answer = raw_input("> ")
        else:
            answer = auto_items[auto_ind]
            auto_ind += 1
        if answer == '0':
            preliminary_flagging(vis=vis, my_line_spws=my_line_spws, 
                                 my_cont_spws=my_cont_spws,
                                 shadow_tolerance=shadow_tolerance,
                                 quackinterval=quackinterval,
                                 config=config)
        elif answer == '1':
            auto_flag_calibrators(vis=vis, primary_cals=primary_cals, 
                                  secondary_cals=secondary_cals)            
        elif answer == '2':
            gen_calibrator_plots(vis=vis, primary_cals=primary_cals, 
                                 secondary_cals=secondary_cals, 
                                 config=config)
        elif answer == '3':
            manual_flag_calibrators(vis=vis, primary_cals=primary_cals, 
                                    secondary_cals=secondary_cals, 
                                    science_targets=science_targets, 
                                    config=config)
        elif answer == '4':
            caltables, calfields = \
              calibrate_calibrators(vis=vis, primary_cals=primary_cals, 
                                    secondary_cals=secondary_cals, 
                                    flux_cals=flux_cals, 
                                    my_line_spws=my_line_spws, 
                                    my_cont_spws=my_cont_spws, 
                                    refant=refant, calwt=calwt,
                                    antpos=antpos, gaincurve=gaincurve,
                                    opacity=opacity, config=config)
        elif answer == '5':
            calibrate_sciencetargets(caltables, calfields, vis=vis,
                                     science_targets=science_targets, calwt=calwt)
        elif answer == '6':
            auto_flag_sciencetargets(vis=vis, science_targets=science_targets)
        elif answer == '7':
            gen_sciencetarget_plots(vis=vis, science_targets=science_targets, 
                                    config=config)
        elif answer == '8':
            manual_flag_sciencetargets(vis=vis, science_targets=science_targets, 
                                       config=config)
        elif answer == '9':
            split_fields(vis=vis, primary_cals=primary_cals, secondary_cals=secondary_cals, 
                         science_targets=science_targets)
        elif answer.lower() == 'q' or answer.lower() == 'quit':
            break
        else:
            print("Input not recognized.")
        if auto_ind >= len(auto_items):
            break
