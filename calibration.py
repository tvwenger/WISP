"""
calibration.py - WISP Calibration Pipeline

Calibrate a measurement set by scripting CASA tasks and generating
diagnostic plots.

Copyright(C) 2018-2019 by
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

Trey V. Wenger August 2019 - V2.0
    Added polarization calibration.
    Re-designed to OOP framework.
"""

import os
import glob
import re
import time
import pickle
import gc

import ConfigParser
import logging
import logging.config

import numpy as np
from scipy.interpolate import interp1d

import __main__ as casa
from recipes.atcapolhelpers import qufromgain

__version__ = '2.0'

# load logging configuration file
logging.config.fileConfig('logging.conf')

# catch re-naming of raw_input to input in Python3
try:
    input = raw_input
except NameError:
    pass

def natural_sort(mylist):
    """
    Natural sort an alphanumeric list

    Inputs: mylist
      mylist :: list of strings
        The list of strings to be sorted

    Returns: sorted_list
      sorted_list :: list of strings
        The sorted list
    """
    # Convert text to integers or lowercase
    convert = \
        lambda text: int(text) if text.isdigit() else text.lower()
    # define the sorting algorithm
    alphanum_key = lambda key: [convert(c) for c in
                                re.split('([0-9]+)', key)]
    # return the sorted list
    return sorted(mylist, key=alphanum_key)

class Calibration:
    """
    The Calibration object handles the calibration steps for a
    measurement set.
    """

    def __init__(self, vis, logger, config,
                 shadow_tolerance=0.0, quack_interval=10.0,
                 antpos=True, gaincurve=True, opacity=True,
                 calpol=False, calwt=True):
        """
        Create a new Calibration object. Get reference antenna,
        create listobs file, find line and continuum spectral windows,
        generate lists of calibrators, and create needed directories.

        Inputs:
          vis :: string
            The masurement set
          logger :: logging.Logger object
            The logging object we're using
          config :: config.ConfigParser object
            The config parser we're using
          shadow_tolerance :: scalar
            The overlap tolerance used for shadow flagging. Flag
            any data with projected baseline separation less than
            r_1 + r_2 - shadow_tolerance
            where r_1 and r_2 are the radii of the antennas.
          quack_interval :: scalar
            The amount of time in seconds to flag at the beginning
            of each scan.
          antpos :: boolean
            if True, compute antenna position corrections
            (only for VLA)
          gaincurve :: boolean
            if True, compute gain curve and antenna efficiency
            corrections (only for VLA)
          opacity :: boolean
            if True, compute opacity corrections
          calpol :: boolean
            if True, calibrate polarization
          calwt :: boolean
            if True, apply calibration weights to data

        Returns: calibration
          calibration :: calibration.Calibration object
            a new Calibration object
        """
        self.vis = vis
        self.logger = logger
        self.config = config
        self.shadow_tolerance = shadow_tolerance
        self.quack_interval = quack_interval
        self.antpos = antpos
        self.gaincurve = gaincurve
        self.opacity = opacity
        self.calpol = calpol
        self.calwt = calwt
        self.logger.info('Initializing Calibration object.')
        #
        # Initialize calibration tables
        #
        # Storage for optional tables (antpos, gaincurve, opacity)
        self.gaintables = []
        self.gainfields = []
        if self.antpos and os.path.isdir('antpos.cal'):
            self.gaintables += ['antpos.cal']
            self.gainfields += ['']
            self.logger.info('Loaded antpos.cal')
        if self.gaincurve and os.path.isdir('gaincurve.cal'):
            self.gaintables += ['gaincurve.cal']
            self.gainfields += ['']
            self.logger.info('Loaded gaincurve.cal')
        if self.opacity and os.path.isdir('opacity.cal'):
            self.gaintables += ['opacity.cal']
            self.gainfields += ['']
            self.logger.info('Loaded opacity.cal')
        # Storage for most recent delays.Kcal table
        self.delays = None
        if os.path.isdir('delays.Kcal'):
            self.delays = 'delays.Kcal'
            self.logger.info('Loaded delays.Kcal')
        # Storage for most recent bandpass.Bcal table
        self.bandpass = None
        check = glob.glob('bandpass.Bcal*')
        if check:
            check.sort()
            self.bandpass = check[-1]
            self.logger.info('Loaded %s', self.bandpass)
        # Storage for most recent phase_int.Gcal table
        self.phase_int = None
        check = glob.glob('phase_int.Gcal*')
        if check:
            check.sort()
            self.phase_int = check[-1]
            self.logger.info('Loaded %s', self.phase_int)
        # Storage for most recent phase_scan.Gcal table
        self.phase_scan = None
        check = glob.glob('phase_scan.Gcal*')
        if check:
            check.sort()
            self.phase_scan = check[-1]
            self.logger.info('Loaded %s', self.phase_scan)
        # Storage for most recent apcal_scan.Gcal table
        self.apcal_scan = None
        check = glob.glob('apcal_scan.Gcal*')
        if check:
            check.sort()
            self.apcal_scan = check[-1]
            self.logger.info('Loaded %s', self.apcal_scan)
        # Storage for flux.cal table
        self.flux = None
        if os.path.isdir('flux.cal'):
            self.flux = 'flux.cal'
            self.logger.info('Loaded flux.cal')
        # Storage for most recent polcal_scan.Dcal table
        self.polcal_scan = None
        check = glob.glob('polcal_scan.Dcal*')
        if check:
            check.sort()
            self.polcal_scan = check[-1]
            self.logger.info('Loaded %s', self.polcal_scan)
        #
        # Get reference antenna
        #
        self.logger.info('Looking for good reference antenna...')
        self.refant = input('Reference Antenna? ')
        if self.refant is None:
            self.logger.critical('No good referance antenna found!')
            raise ValueError('No good referance antenna found!')
        self.logger.info('Done. Found reference antenna: %s',
                         self.refant)
        #
        # Generate listobs file
        #
        listfile = 'listobs.txt'
        if not os.path.isfile(listfile):
            self.logger.info('Generating listobs file...')
            casa.listobs(vis=self.vis, listfile=listfile)
            self.logger.info('Done.')
        #
        # Get continuum and line spws from configuration file
        #
        self.cont_spws = self.config.get('Spectral Windows',
                                         'Continuum')
        self.line_spws = self.config.get('Spectral Windows', 'Line')
        self.logger.info('Found continuum spws: %s', self.cont_spws)
        self.logger.info('Found line spws: %s', self.line_spws)
        #
        # get field names
        #
        self.logger.info('Looking for field names...')
        self.all_fields = casa.vishead(
            vis=self.vis, mode='get', hdkey='field')[0]
        self.all_fields = list(set(self.all_fields))
        self.logger.info('Found fields: %s',
                         ', '.join(self.all_fields))
        #
        # Get primary calibrator fields if they are not in config
        #
        config_pri_cals = self.config.get('Calibrators',
                                          'Primary Calibrators')
        if config_pri_cals == '':
            self.pri_cals = []
            self.logger.info('Looking for primary calibrators...')
            with open(listfile, 'r') as fin:
                for line in fin:
                    for field in self.all_fields:
                        if ((field in line) and
                                ('CALIBRATE_BANDPASS' in line) and
                                (field not in self.pri_cals)):
                            self.pri_cals.append(field)
            self.logger.info('Done')
        else:
            self.pri_cals = [field for field in
                             config_pri_cals.splitlines()
                             if field in self.all_fields]
        self.logger.info('Primary calibrators: %s',
                         ', '.join(self.pri_cals))
        #
        # Get Secondary calibrator fields if they are not in config
        #
        config_sec_cals = self.config.get('Calibrators',
                                          'Secondary Calibrators')
        if config_sec_cals == '':
            self.sec_cals = []
            self.logger.info('Looking for secondary calibrators...')
            with open(listfile, 'r') as fin:
                for line in fin:
                    for field in self.all_fields:
                        if ((field in line) and
                                (('CALIBRATE_AMPLI' in line) or
                                 ('CALIBRATE_PHASE' in line)) and
                                (field not in self.sec_cals)):
                            self.sec_cals.append(field)
            self.logger.info('Done')
        else:
            self.sec_cals = [field for field in
                             config_sec_cals.splitlines()
                             if field in self.all_fields]
        self.logger.info('Secondary calibrators: %s',
                         ', '.join(self.sec_cals))
        #
        # Get flux calibrator fields if they are not in config
        #
        config_flux_cals = self.config.get('Calibrators',
                                           'Flux Calibrators')
        if config_flux_cals == '':
            self.flux_cals = []
            self.logger.info('Looking for flux calibrators...')
            with open(listfile, 'r') as fin:
                for line in fin:
                    for field in self.all_fields:
                        if ((field in line) and
                                ('CALIBRATE_FLUX' in line) and
                                (field not in self.flux_cals)):
                            self.flux_cals.append(field)
            self.logger.info('Done')
        else:
            self.flux_cals = [field for field in
                              config_flux_cals.splitlines()
                              if field in self.all_fields]
        self.logger.info('Flux calibrators: %s',
                         ', '.join(self.flux_cals))
        #
        # Check that flux calibrators are in primary calibrator list
        # if not, add them
        #
        for flux_cal in self.flux_cals:
            if flux_cal not in self.pri_cals:
                self.pri_cals.append(flux_cal)
        #
        # Get polarization calibrators from config
        #
        config_pol_cals = self.config.get('Calibrators',
                                          'Polarization Calibrators')
        self.pol_cals = [field for field in
                         config_pol_cals.splitlines()
                         if field in self.all_fields]
        #
        # check that polarization calibrators are in primary list
        #
        for pol_cal in self.pol_cals:
            if pol_cal not in self.pri_cals:
                raise ValueError('Polarization Calibrators must be a '
                                 'Primary Calibrator')
        self.logger.info('Polarization calibrators: %s',
                         ', '.join(self.pol_cals))
        #
        # Get science targets
        #
        self.sci_targets = []
        self.logger.info('Looking for science targets...')
        for field in self.all_fields:
            if field not in self.pri_cals+self.sec_cals:
                self.sci_targets.append(field)
        self.logger.info('Done')
        self.logger.info('Science targets: %s',
                         ', '.join(self.sci_targets))
        #
        # Determine if we are doing parallactic angle correction
        #
        self.parang = self.config.getboolean(
            'Polarization', 'Parallactic Angle Correction')
        #
        # create directories for figures
        #
        if not os.path.isdir('calibrator_figures'):
            os.makedirs('calibrator_figures')
        if not os.path.isdir('science_figures'):
            os.makedirs('science_figures')
        if not os.path.isdir('plotcal_figures'):
            os.makedirs('plotcal_figures')

    def save_flags(self, label):
        """
        Save the flag state with a given label and the current
        time.

        Inputs:
          label :: string
            A label to add to the version name

        Returns: Nothing
        """
        self.logger.info('Saving flag state...')
        cur_time = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
        versionname = '{0} {1}'.format(label, cur_time)
        casa.flagmanager(vis=self.vis, mode='save',
                         versionname=versionname)
        self.logger.info('Done')

    def interpolate_channels(self, spws, chans):
        """
        Edit the measurement set to replace bad channels with
        interpolated values. Simple linear interpolation between
        neighbors, in both phase and amplitude.

        Inputs:
          spws :: list of integers
            The spectral windows to use
          chans :: list of integers
            The channels through which to interpolate

        Returns: Nothing
        """
        bad_chans = np.array(chans)
        for spw in spws:
            self.logger.info('Working on spw %d', spw)
            casa.ms.open(self.vis, nomodify=False)
            casa.ms.selectinit(datadescid=spw)
            self.logger.info('Reading data...')
            data = casa.ms.getdata(['data'])
            self.logger.info('Done.')
            chans = np.arange(data['data'].shape[1])
            mask = np.zeros(data['data'].shape[1], dtype=bool)
            mask[bad_chans] = True
            #
            # Interpolate amplitude
            #
            self.logger.info('Interpolating amplitude...')
            amp = np.abs(data['data'])
            amp_interp = interp1d(chans[~mask], amp[:, ~mask, :],
                                  axis=1)
            amp[:, mask, :] = amp_interp(chans[mask])
            self.logger.info('Done.')
            #
            # Interpolate phase
            #
            self.logger.info('Interpolating phase...')
            phase = np.unwrap(np.angle(data['data']))
            phase_interp = interp1d(chans[~mask], phase[:, ~mask, :],
                                    axis=1)
            phase[:, mask, :] = phase_interp(chans[mask])
            self.logger.info('Done.')
            #
            # Save
            #
            self.logger.info('Updating visibilities...')
            data['data'] = amp*np.cos(phase) + 1.j*amp*np.sin(phase)
            self.logger.info('Done.')
            self.logger('Saving to measurement set...')
            casa.ms.putdata(data)
            casa.ms.done()
            self.logger.info('Done.')
            del data
            del amp
            del amp_interp
            del phase
            del phase_interp
            gc.collect()

    def preliminary_flagging(self):
        """
        Perform preliminary flagging: shadowed antennas, quack,
        flags from configuration file, interpolations from
        configuration file, then extend flags as necessary.

        Inputs: Nothing

        Returns: Nothing
        """
        self.save_flags('starting_flags')
        #
        # Flag shadowed antennas
        #
        self.logger.info('Flagging shadowed antennas...')
        casa.flagdata(vis=self.vis, mode='shadow',
                      tolerance=self.shadow_tolerance,
                      flagbackup=False, extendflags=False)
        self.logger.info('Done.')
        #
        # Flag the beginning of each scan
        #
        self.logger.info('Flagging the beginning of each scan (quack)...')
        casa.flagdata(vis=self.vis, mode='quack',
                      quackinterval=self.quack_interval,
                      flagbackup=False, extendflags=False)
        self.logger.info('Done.')
        #
        # Flag scans from configuration file
        #
        scan = self.config.get('Flags', 'Scan')
        if scan != '':
            self.logger.info('Flagging scans from configuration file: '
                             '%s', scan)
            casa.flagdata(vis=self.vis, mode='manual', scan=scan,
                          flagbackup=False, extendflags=False)
            self.logger.info('Done.')
        #
        # Flag antennas from configuration file
        #
        antenna = self.config.get('Flags', 'Antenna')
        if antenna != '':
            self.logger.info('Flagging antennas from configuration '
                             'file: %s', antenna)
            casa.flagdata(vis=self.vis, mode='manual', antenna=antenna,
                          flagbackup=False, extendflags=False)
            self.logger.info('Done.')
        #
        # Flag spectral windows from configuration file
        #
        spws = self.config.get('Flags', 'Spectral Window')
        if spws != '':
            self.logger.info('Flagging spectral windows from '
                             'configuration file: %s', spws)
            casa.flagdata(vis=self.vis, mode='manual', spw=spws,
                          flagbackup=False, extendflags=False)
            self.logger.info('Done.')
        #
        # Flag line channels from configuration file
        #
        badchans = self.config.get('Flags', 'Line Channels')
        badchans = badchans.split(',')
        if badchans[0] != '':
            self.logger.info('Flagging line channels from configuration '
                             'file: %s', ';'.join(badchans))
            badchans = ';'.join(badchans)
            line_spws = ','.join([i+':'+badchans for i in
                                  self.line_spws.split(',')])
            casa.flagdata(vis=self.vis, mode='manual', spw=line_spws,
                          flagbackup=False, extendflags=False)
            self.logger.info('Done.')
        #
        # Flag continuum channels from configuration file
        #
        badchans = self.config.get('Flags', 'Continuum Channels')
        badchans = badchans.split(',')
        if badchans[0] != '':
            self.logger.info('Flagging continuum channels from '
                             'configuration file: %s',
                             ';'.join(badchans))
            badchans = ';'.join(badchans)
            cont_spws = ','.join([i+':'+badchans for i in
                                  self.cont_spws.split(',')])
            casa.flagdata(vis=self.vis, mode='manual', spw=cont_spws,
                          flagbackup=False, extendflags=False)
            self.logger.info('Done.')
        #
        # Interpolate through bad line channels
        #
        badchans = self.config.get('Interpolate', 'Line Channels')
        badchans = badchans.split(',')
        if badchans[0] != '':
            self.logger.info('Interpolating through line channels from '
                             'configuration file: %s',
                             ';'.join(badchans))
            line_spws = [int(i) for i in self.line_spws.split(',')]
            badchans = np.array([int(i) for i in badchans])
            self.interpolate_channels(line_spws, badchans)
            self.logger.info('Done.')
        #
        # Interpolate through bad continuum channels
        #
        badchans = self.config.get('Interpolate', 'Continuum Channels')
        badchans = badchans.split(',')
        if badchans[0] != '':
            self.logger.info('Interpolating through continuum channels '
                             'from configuration file: %s',
                             ';'.join(badchans))
            cont_spws = [int(i) for i in self.cont_spws.split(',')]
            badchans = np.array([int(i) for i in badchans])
            self.interpolate_channels(cont_spws, badchans)
            self.logger.info('Done')
        #
        # Extend the flags
        #
        self.logger.info('Extending flags...')
        casa.flagdata(vis=self.vis, mode='extend', extendpols=True,
                      growtime=75.0, growfreq=50.0, growaround=True,
                      flagbackup=False)
        self.logger.info('Done.')
        self.save_flags('preliminary')

    def auto_flag_calibrators(self):
        """
        Perform automatic flagging of calibrators using rflag on
        calibrated data, then extend the flags.

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # check if calibrators have corrected datacolumn
        #
        field = ','.join(self.pri_cals+self.sec_cals)
        stat = None
        self.logger.info('Checking if ms contains corrected data '
                         'column...')
        stat = casa.visstat(vis=self.vis, field=field, spw='0',
                            datacolumn='corrected')
        if stat is None:
            self.logger.critical('Done. ms does not contain corrected '
                                 'data column. Skipping.')
            return
        self.logger.info('Done. ms does contain corrected data column.')
        #
        # Run rflag on calibrated data
        #
        self.logger.info('Running rflag on corrected data column...')
        casa.flagdata(vis=self.vis, mode='rflag', field=field,
                      flagbackup=False, datacolumn='corrected',
                      extendflags=False)
        self.logger.info('Done.')
        #
        # Extend the flags and save
        #
        self.logger.info('Extending flags...')
        casa.flagdata(vis=self.vis, mode='extend', extendpols=True,
                      growtime=75.0, growfreq=50.0, growaround=True,
                      flagbackup=False)
        self.logger.info('Done.')
        self.save_flags('autoflag')

    def calibrator_plots(self):
        """
        Generate diagnostic visiblity plots for calibrators.

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # check if calibrators have corrected datacolumn
        #
        field = ','.join(self.pri_cals+self.sec_cals)
        stat = None
        self.logger.info('Checking if ms contains corrected data '
                         'column...')
        stat = casa.visstat(vis=self.vis, spw='0', field=field,
                            datacolumn='corrected')
        if stat is None:
            self.logger.info('Done. ms does not contain corrected '
                             'data column.')
            datacolumn = 'data'
        else:
            self.logger.info('Done. ms does contain corrected data '
                             'column.')
            datacolumn = 'corrected'
        #
        # Generate the plots
        #
        self.logger.info('Generating plots for manual inspection...')
        corr = self.config.get('Polarization', 'Polarization')
        plotnum = 0
        plots = []
        for field in self.pri_cals+self.sec_cals:
            #
            # Phase vs. Amplitude
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'calibrator_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='amp', yaxis='phase',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=corr, title=title,
                        plotfile=plotfile, overwrite=True,
                        showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'amp', 'yaxis':'phase',
                          'avgtime':'', 'avgchannel':''})
            plotnum += 1
            #
            # Amplitude vs UV-distance (in wavelength units)
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'calibrator_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='uvwave', yaxis='amp',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=corr, title=title,
                        plotfile=plotfile, overwrite=True,
                        showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'uvwave', 'yaxis':'amp',
                          'avgtime':'', 'avgchannel':''})
            plotnum += 1
            #
            # Amplitude vs Time
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'calibrator_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='time', yaxis='amp',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=corr, avgchannel='1e7',
                        title=title, plotfile=plotfile,
                        overwrite=True, showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'time', 'yaxis':'amp',
                          'avgtime':'', 'avgchannel':'1e7'})
            plotnum += 1
            #
            # Phase vs Time
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'calibrator_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='time', yaxis='phase',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=corr, avgchannel='1e7',
                        title=title, plotfile=plotfile,
                        overwrite=True, showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'time', 'yaxis':'phase',
                          'avgtime':'', 'avgchannel':'1e7'})
            plotnum += 1
            #
            # Amplitude vs Channel
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'calibrator_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='channel', yaxis='amp',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=corr, avgtime='1e7',
                        title=title, plotfile=plotfile,
                        overwrite=True, showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'channel', 'yaxis':'amp',
                          'avgtime':'1e7', 'avgchannel':''})
            plotnum += 1
            #
            # Phase vs Channel
            # 
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'calibrator_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='channel', yaxis='phase',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=corr, avgtime='1e7',
                        title=title, plotfile=plotfile,
                        overwrite=True, showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'channel', 'yaxis':'phase',
                          'avgtime':'1e7', 'avgchannel':''})
            plotnum += 1
        self.logger.info('Done.')
        #
        # Generate PDF to display plots
        #
        self.logger.info('Generating PDF...')
        num_plots = plotnum
        iplot = 0
        with open('calibrator_figures.tex', 'w') as fout:
            fout.write(r'\documentclass{article}'+'\n')
            fout.write(r'\usepackage{graphicx}'+'\n')
            fout.write(r'\usepackage[margin=0.1cm]{geometry}'+'\n')
            fout.write(r'\begin{document}'+'\n')
            fout.write(r'\begin{figure}'+'\n')
            fout.write(r'\centering'+'\n')
            for plotnum in range(num_plots):
                fnames = glob.glob(
                    'calibrator_figures/{0}_*.png'.format(plotnum))
                fnames = natural_sort(fnames)
                for fname in fnames:
                    if iplot > 0 and iplot % 6 == 0:
                        fout.write(r'\end{figure}'+'\n')
                        fout.write(r'\clearpage'+'\n')
                        fout.write(r'\begin{figure}'+'\n')
                        fout.write(r'\centering'+'\n')
                    elif iplot > 0 and iplot % 2 == 0:
                        fout.write(r'\end{figure}'+'\n')
                        fout.write(r'\begin{figure}'+'\n')
                        fout.write(r'\centering'+'\n')
                    fout.write(r'\includegraphics[width=0.45\textwidth]'
                               '{'+fname+'}\n')
                    iplot += 1
            fout.write(r'\end{figure}'+'\n')
            fout.write(r'\end{document}'+'\n')
        os.system('pdflatex -interaction=batchmode calibrator_figures.tex')
        self.logger.info('Done.')
        #
        # Save plot list to a pickle object
        #
        self.logger.info('Saving plot list to pickle...')
        with open('calibrator_plots.pkl', 'w') as fout:
            pickle.dump(plots, fout)
        self.logger.info('Done.')

    def flag(self, fields):
        """
        Interactively flag some data.

        Inputs:
          fields :: list of strings
            The fields to flag

        Returns: Nothing
        """
        #
        # Build list of flag commands
        #
        flag_commands = []
        while True:
            #
            # Prompt user for attributes to flag
            #
            print('Field? Empty = {0}'.format(','.join(fields)))
            field = input()
            if field == '':
                field = ','.join(fields)
            if field not in ','.join(fields):
                print('{0} is not in {1}'.format(field, fields))
                continue
            print('Scan? Empty = all scans (ex. 0 to flag scan 0, 1~3 '
                  'to flag scans 1, 2, and 3)')
            scan = input()
            print('Spectral window and channels? Empty = all spws '
                  '(ex. 2:100~150;200:250,5:1200~1210 is spw 2, chans '
                  '100 to 150 and 200 to 250, and spw 5 chans 1200 to '
                  '1210)')
            spw = input()
            print('Time range? Empty = all times (ex. '
                  '10:23:45~10:23:55)')
            timerange = input()
            print('Antenna or baseline? Empty = all antennas/baselines '
                  '(ex. CA01 to flag ant 1 or CA01&CA02 to flag 1-2 '
                  'baseline)')
            antenna = input()
            print('Correlation? Empty = all correlations (i.e. XX,YY)')
            correlation = input()
            #
            # Build flag command
            #
            parts = []
            if field:
                parts.append("field='{0}'".format(field))
            if scan:
                parts.append("scan='{0}'".format(scan))
            if spw:
                parts.append("spw='{0}'".format(spw))
            if timerange:
                parts.append("timerange='{0}'".format(timerange))
            if antenna:
                parts.append("antenna='{0}'".format(antenna))
            if correlation:
                parts.append("correlation='{0}'".format(correlation))
            flag_commands.append(' '.join(parts))
            #
            # Confirm with user, or append more flag commands
            #
            print('Will execute:')
            print("flagdata(vis='{0}',mode='list',flagbackup=False,"
                  "extendflags=False,".format(self.vis))
            for icmd, cmd in enumerate(flag_commands):
                if len(flag_commands) == 1:
                    print("         inpfile=[\"{0}\"])".format(cmd))
                elif icmd == 0:
                    print("         inpfile=[\"{0}\",".format(cmd))
                elif icmd == len(flag_commands)-1:
                    print("                  \"{0}\"])".format(cmd))
                else:
                    print("                  \"{0}\",".format(cmd))
            print('Proceed [y/n] or add another flag command [a]?')
            answer = input()
            #
            # Execute flag command
            #
            if answer.lower() == 'y':
                self.logger.info("Executing:")
                self.logger.info("flagdata(vis='%s',mode='list',"
                                 "flagbackup=False,extendflags=False,",
                                 self.vis)
                for icmd, cmd in enumerate(flag_commands):
                    if len(flag_commands) == 1:
                        self.logger.info("         inpfile=[\"%s\"])",
                                         cmd)
                    elif icmd == 0:
                        self.logger.info("         inpfile=[\"%s\",",
                                         cmd)
                    elif icmd == len(flag_commands)-1:
                        self.logger.info("                  \"%s\"])",
                                         cmd)
                    else:
                        self.logger.info("                  \"%s\",",
                                         cmd)
                #
                # Save flag command to manual flags list
                #
                with open('manual_flags.txt', 'a') as fout:
                    cur_time = time.strftime('%Y-%m-%d %H:%M:%S',
                                             time.gmtime())
                    for cmd in flag_commands:
                        fout.write('{0}: {1}'.format(cur_time, cmd)+
                                   '\n')
                #
                # Execute
                #
                casa.flagdata(vis=self.vis, mode='list',
                              flagbackup=False, extendflags=False,
                              inpfile=flag_commands)
                break
            #
            # Append another flag command
            #
            elif answer.lower() == 'a':
                continue
            #
            # Quit flagging
            #
            else:
                print('Aborting...')
                break

    def manual_flag_calibrators(self):
        """
        Interactively plot and flag the calibrators.

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # check if calibrators have corrected datacolumn
        #
        field = ','.join(self.pri_cals+self.sec_cals)
        stat = None
        self.logger.info('Checking if ms contains corrected data '
                         'column...')
        stat = casa.visstat(vis=self.vis, spw='0', field=field,
                            datacolumn='corrected')
        if stat is None:
            self.logger.info('Done. ms does not contain corrected '
                             'data column.')
            datacolumn = 'data'
        else:
            self.logger.info('Done. ms does contain corrected data '
                             'column.')
            datacolumn = 'corrected'
        #
        # Read the plot list from the pickle object
        #
        self.logger.info('Reading plot list from pickle...')
        with open('calibrator_plots.pkl', 'r') as fin:
            plots = pickle.load(fin)
        num_plots = len(plots)
        self.logger.info('Done.')
        #
        # Display menu option to user
        #
        corr = self.config.get('Polarization', 'Polarization')
        self.logger.info('Please inspect calibrator_plots.pdf then '
                         'perform manual flagging.')
        while True:
            print('f - flag some data')
            print('plot id number - generate interactive version of plot with this id')
            print('quit - end this flagging session')
            answer = input()
            #
            # Flag some data
            #
            if answer.lower() == 'f':
                # if we are going to flag something in all fields for
                # the calibrators, we might as well flag it in the
                # science targets too
                self.flag(self.all_fields)
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
                    print('Invalid Plot ID')
                    continue
                if plotid >= num_plots:
                    print('Invalid Plot ID')
                    continue
                title = 'PlotID: {0} Field: {1}'.\
                    format(plotid, plots[plotid]['field'])
                casa.plotms(vis=self.vis,
                            xaxis=plots[plotid]['xaxis'],
                            yaxis=plots[plotid]['yaxis'],
                            field=plots[plotid]['field'],
                            ydatacolumn=datacolumn,
                            iteraxis='spw', coloraxis='baseline',
                            correlation=corr, title=title,
                            avgchannel=plots[plotid]['avgchannel'],
                            avgtime=plots[plotid]['avgtime'])
        #
        # Save the flags
        #
        self.save_flags('manualflag')

    def calibration_tables(self):
        """
        Generate gain, bandpass, and polarization leakage
        calibration tables.

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # bandpass calibration for continuum spws. Combine all scans,
        # average some channels as defined in configuration file
        #
        caltable = 'bandpass.Bcal0'
        gaintables = self.gaintables + [self.delays, self.phase_int]
        self.logger.info('Calculating bandpass calibration table for '
                         'primary calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        chan_avg = self.config.get('Bandpass Channel Average',
                                   'Continuum Channels')
        if chan_avg == '':
            solint = 'inf'
        else:
            solint = 'inf,{0}chan'.format(chan_avg)
        field = ','.join(self.pri_cals)
        casa.bandpass(vis=self.vis, caltable=caltable, field=field,
                      spw=self.cont_spws, refant=self.refant,
                      solint=solint, combine='scan', solnorm=True,
                      minblperant=1, parang=self.parang,
                      gaintable=gaintables)
        #
        # bandpass calibration for line spws. Combine all scans,
        # average some channels as defined in configuration file,
        # append to continuum channel bandpass calibration table
        #
        chan_avg = self.config.get('Bandpass Channel Average',
                                   'Line Channels')
        if chan_avg == '':
            solint = 'inf'
        else:
            solint = 'inf,{0}chan'.format(chan_avg)
        casa.bandpass(vis=self.vis, caltable=caltable, field=field,
                      spw=self.line_spws, refant=self.refant,
                      solint=solint, combine='scan', solnorm=True,
                      minblperant=1, append=True, parang=self.parang,
                      gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with bandpass calibration')
            raise ValueError('Problem with bandpass calibration!')
        self.bandpass = caltable
        self.logger.info('Done.')
        #
        # integration timescale phase corrections for all calibrators
        # required for accurate amplitude calibration
        #
        caltable = 'phase_int.Gcal1'
        gaintables = self.gaintables + [self.delays, self.bandpass]
        self.logger.info('Re-calculating the phase calibration table '
                         'on integration timescales for all '
                         'calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     solint='int', calmode='p', refant=self.refant,
                     gaintype='G', minsnr=2.0, minblperant=1,
                     parang=self.parang, gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with integration-timescale '
                                 'phase calibration')
            raise ValueError('Problem with integration-timescale '
                             'phase calibration!')
        self.phase_int = caltable
        self.logger.info('Done.')
        #
        # scan timescale phase corrections for all calibrators
        # required to apply to science targets
        #
        caltable = 'phase_scan.Gcal0'
        gaintables = self.gaintables + [self.delays, self.bandpass]
        self.logger.info('Calculating the phase calibration table on '
                         'scan timescales for all calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     solint='inf', calmode='p', refant=self.refant,
                     gaintype='G', minsnr=2.0, minblperant=1,
                     parang=self.parang, gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with scan-timescale phase '
                                 'calibration')
            raise ValueError('Problem with scan-timescale phase '
                             'calibration!')
        self.phase_scan = caltable
        self.logger.info('Done.')
        #
        # scan timescale amplitude corrections using
        # integration timescale phase calibration
        #
        caltable = 'apcal_scan.Gcal0'
        gaintables = self.gaintables + [
            self.delays, self.bandpass, self.phase_int]
        self.logger.info('Calculating the amplitude calibration table '
                         'on scan timescales for all calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     solint='inf', calmode='ap', refant=self.refant,
                     minsnr=2.0, minblperant=1, parang=self.parang,
                     gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with amplitude calibration')
            raise ValueError('Problem with amplitude calibration!')
        self.apcal_scan = caltable
        self.logger.info('Done.')
        #
        # Done if not doing polarization calibration
        #
        if not self.calpol:
            return
        #
        # Polarization leakage calibration on scan timescales
        # Only continuum spws
        #
        caltable = 'polcal_scan.Dcal0'
        gaintables = self.gaintables + [
            self.delays, self.bandpass, self.phase_int,
            self.apcal_scan]
        self.logger.info('Calculating the polarization leakage '
                         'calibration table on scan timescales for '
                         'polarization calibrators in '
                         'continuum spectral windows...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.pol_cals)
        casa.polcal(vis=self.vis, caltable=caltable, field=field,
                    spw=self.cont_spws, solint='inf', poltype='D',
                    refant=self.refant, minsnr=2.0, minblperant=1,
                    gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with polarization '
                                 'leakage calibration')
            raise ValueError('Problem with polarization leakage '
                             'calibration')
        self.polcal_scan = caltable
        self.logger.info('Done.')

    def calibration_tables_post_polcal(self):
        """
        Re-generate gain, bandpass, and polarization leakage
        calibration tables, now with initial estimates for source
        polarizations.

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # Get source polarization estimates
        #
        smodels = {}
        # need to get the field ID numbers, so get all fields in
        # order
        field_order = casa.vishead(vis=self.vis, mode='get',
                                   hdkey='field')[0]
        for field in self.pri_cals+self.sec_cals:
            # skip if this is a polarization calibrator
            if field in self.pol_cals:
                continue
            # get fieldid(s)
            fieldids = [i for i in range(len(field_order))
                        if field_order[i] == field]
            # get QU first-order correction.
            qu = qufromgain(self.apcal_scan, fieldids=fieldids)
            # average all spectral windows.
            estimate_q = np.mean([qu[fieldid][0] for fieldid in
                                  fieldids])
            estimate_u = np.mean([qu[fieldid][1] for fieldid in
                                  fieldids])
            # set polarization solution
            smodels[field] = [1, estimate_q, estimate_u, 0]
        #
        # bandpass calibration for continuum spws. Combine all scans,
        # average some channels as defined in configuration file,
        # and apply polarization leakage calibration
        #
        caltable = 'bandpass.Bcal1'
        gaintables = self.gaintables + [
            self.delays, self.phase_int, self.apcal_scan,
            self.polcal_scan]
        self.logger.info('Calculating bandpass calibration table for '
                         'primary calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        chan_avg = self.config.get('Bandpass Channel Average',
                                   'Continuum Channels')
        if chan_avg == '':
            solint = 'inf'
        else:
            solint = 'inf,{0}chan'.format(chan_avg)
        #
        # First the polarization calibrator
        #
        field = ','.join(self.pol_cals)
        casa.bandpass(vis=self.vis, caltable=caltable, field=field,
                      spw=self.cont_spws, refant=self.refant,
                      solint=solint, combine='scan', solnorm=True,
                      minblperant=1, parang=self.parang,
                      gaintable=gaintables)
        #
        # Then the other primary calibrators
        #
        for field in self.pri_cals:
            if field in self.pol_cals:
                continue
            casa.bandpass(vis=self.vis, caltable=caltable,
                          field=field, spw=self.cont_spws,
                          refant=self.refant, solint=solint,
                          combine='scan', solnorm=True, minblperant=1,
                          parang=self.parang, gaintable=gaintables,
                          smodel=smodels[field], append=True)
        #
        # bandpass calibration for line spws. Combine all scans,
        # average some channels as defined in configuration file,
        # append to continuum channel bandpass calibration table.
        # No polarization leakage calibration
        #
        gaintables = self.gaintables + [
            self.delays, self.phase_int, self.apcal_scan]
        chan_avg = self.config.get('Bandpass Channel Average',
                                   'Line Channels')
        if chan_avg == '':
            solint = 'inf'
        else:
            solint = 'inf,{0}chan'.format(chan_avg)
        field = ','.join(self.pri_cals)
        casa.bandpass(vis=self.vis, caltable=caltable, field=field,
                      spw=self.line_spws, refant=self.refant,
                      solint=solint, combine='scan', solnorm=True,
                      minblperant=1, append=True, parang=self.parang,
                      gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with bandpass calibration')
            raise ValueError('Problem with bandpass calibration!')
        self.bandpass = caltable
        self.logger.info('Done.')
        #
        # integration timescale phase corrections for all calibrators
        # required for accurate amplitude calibration. First for the
        # continuum spectral windows including polarization leakage
        # calibration
        #
        caltable = 'phase_int.Gcal2'
        gaintables = self.gaintables + [
            self.delays, self.polcal_scan, self.bandpass]
        self.logger.info('Re-calculating the phase calibration table '
                         'on integration timescales for all '
                         'calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        #
        # First the polarization calibrator
        #
        field = ','.join(self.pol_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     spw=self.cont_spws, solint='int', calmode='p',
                     refant=self.refant, gaintype='G', minsnr=2.0,
                     minblperant=1, parang=self.parang,
                     gaintable=gaintables)
        #
        # Then the others
        #
        for field in self.pri_cals+self.sec_cals:
            if field in self.pol_cals:
                continue
            casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                         spw=self.cont_spws, solint='int',
                         calmode='p', refant=self.refant,
                         gaintype='G', minsnr=2.0, minblperant=1,
                         parang=self.parang, gaintable=gaintables,
                         smodel=smodels[field], append=True)
        #
        # Now the line spectral windows
        #
        gaintables = self.gaintables + [self.delays, self.bandpass]
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     spw=self.line_spws, solint='int', calmode='p',
                     refant=self.refant, gaintype='G', minsnr=2.0,
                     minblperant=1, parang=self.parang,
                     gaintable=gaintables, append=True)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with integration-timescale '
                                 'phase calibration')
            raise ValueError('Problem with integration-timescale '
                             'phase calibration!')
        self.phase_int = caltable
        self.logger.info('Done.')
        #
        # scan timescale phase corrections for all calibrators. First
        # for the continuum spectral windows including polarization
        # leakage calibration
        #
        caltable = 'phase_scan.Gcal1'
        gaintables = self.gaintables + [
            self.delays, self.polcal_scan, self.bandpass]
        self.logger.info('Re-calculating the phase calibration table '
                         'on scan timescales for all '
                         'calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        #
        # First the polarization calibrator
        #
        field = ','.join(self.pol_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     spw=self.cont_spws, solint='inf', calmode='p',
                     refant=self.refant, gaintype='G', minsnr=2.0,
                     minblperant=1, parang=self.parang,
                     gaintable=gaintables)
        #
        # Then the others
        #
        for field in self.pri_cals+self.sec_cals:
            if field in self.pol_cals:
                continue
            casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                         spw=self.cont_spws, solint='inf',
                         calmode='p', refant=self.refant,
                         gaintype='G', minsnr=2.0, minblperant=1,
                         parang=self.parang, gaintable=gaintables,
                         smodel=smodels[field], append=True)
        #
        # Now the line spectral windows
        #
        gaintables = self.gaintables + [self.delays, self.bandpass]
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     spw=self.line_spws, solint='inf', calmode='p',
                     refant=self.refant, gaintype='G', minsnr=2.0,
                     minblperant=1, parang=self.parang,
                     gaintable=gaintables, append=True)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with scan-timescale '
                                 'phase calibration')
            raise ValueError('Problem with scan-timescale '
                             'phase calibration!')
        self.phase_scan = caltable
        self.logger.info('Done.')
        #
        # scan timescale amplitude corrections using
        # integration timescale phase calibration
        #
        caltable = 'apcal_scan.Gcal1'
        gaintables = self.gaintables + [
            self.delays, self.polcal_scan, self.bandpass,
            self.phase_int]
        self.logger.info('Calculating the amplitude calibration table '
                         'on scan timescales for all calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        #
        # First the polarization calibrator
        #
        field = ','.join(self.pol_cals)
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     spw=self.cont_spws, solint='inf', calmode='ap',
                     refant=self.refant, minsnr=2.0, minblperant=1,
                     parang=self.parang, gaintable=gaintables)
        #
        # Now the other calibrators, which we append to the caltable
        #
        for field in self.pri_cals+self.sec_cals:
            if field in self.pol_cals:
                continue
            casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                         spw=self.cont_spws, solint="inf",
                         calmode="ap", refant=self.refant, minsnr=2.0,
                         minblperant=1, parang=self.parang,
                         gaintable=gaintables, smodel=smodels[field],
                         append=True)
        #
        # Now the line spectral windows
        #
        gaintables = self.gaintables + [
            self.delays, self.bandpass, self.phase_int]
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     spw=self.line_spws, solint="inf", calmode="ap",
                     refant=self.refant, minsnr=2.0, minblperant=1,
                     parang=self.parang, gaintable=gaintables,
                     append=True)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with amplitude calibration')
            raise ValueError('Problem with amplitude calibration!')
        self.apcal_scan = caltable
        self.logger.info('Done.')
        #
        # Polarization leakage calibration on scan timescales
        #
        caltable = 'polcal_scan.Dcal1'
        gaintables = self.gaintables + [
            self.delays, self.bandpass, self.phase_int,
            self.apcal_scan]
        self.logger.info('Calculating the polarization leakage '
                         'calibration table on scan timescales for '
                         'polarization calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.pol_cals)
        casa.polcal(vis=self.vis, caltable=caltable, field=field,
                    spw=self.cont_spws, solint='inf', poltype='D',
                    refant=self.refant, minsnr=2.0, minblperant=1,
                    gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with polarization '
                                 'leakage calibration')
            raise ValueError('Problem with polarization leakage '
                             'calibration')
        self.polcal_scan = caltable
        self.logger.info('Done.')

    def calibrate_calibrators(self):
        """
        Calculate calibration solutions and apply solutions to
        the calibrators.

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # Reset calibration tables
        #
        self.gaintables = []
        self.gainfields = []
        self.delays = None
        self.bandpass = None
        self.phase_int = None
        self.phase_scan = None
        self.apcal_scan = None
        self.flux = None
        self.polcal_scan = None
        #
        # set the model for the flux calibrators
        #
        self.logger.info('Setting the flux calibrator models...')
        for flux_cal in self.flux_cals:
            #
            # If flux calibrator model is supplied in config, use that
            #
            manual_flux_cals = self.config.get(
                'Flux Calibrator Models', 'Name')
            manual_flux_cals = manual_flux_cals.splitlines()
            if flux_cal in manual_flux_cals:
                # get index of flux_cal in manual_flux_cals
                flux_idx = manual_flux_cals.index(flux_cal)
                # get reference frequency
                reffreq = self.config.get(
                    'Flux Calibrator Models', 'Reference Frequency')
                reffreq = reffreq.splitlines()[flux_idx]
                # get fluxdensity and convert to proper units
                # Stokes Q, U, and V are assumed 0
                fluxdensity = self.config.get(
                    'Flux Calibrator Models', 'Log Flux Density')
                fluxdensity = fluxdensity.splitlines()
                fluxdensity = [10.**float(fluxdensity[flux_idx]),
                               0., 0., 0.]
                # get spectral index coefficients
                spix = self.config.get(
                    'Flux Calibrator Models',
                    'Spectral Index Coefficients')
                spix = spix.splitlines()[flux_idx]
                spix = [float(i) for i in spix.split(',')]
                # Run setjy in manual mode
                casa.setjy(vis=self.vis, field=flux_cal,
                           scalebychan=True, standard='manual',
                           fluxdensity=fluxdensity, spix=spix,
                           reffreq=reffreq)
            #
            # Otherwise, use CASA model
            #
            else:
                casa.setjy(vis=self.vis, field=flux_cal,
                           scalebychan=True)
        self.logger.info('Done.')
        #
        # Correct antenna positions
        #
        if self.antpos:
            caltable = 'antpos.cal'
            self.logger.info('Calculating antenna position '
                             'correction...')
            if os.path.isdir(caltable):
                casa.rmtables(caltable)
            casa.gencal(vis=self.vis, caltable=caltable,
                        caltype='antpos')
            if os.path.isdir(caltable):
                self.gaintables += [caltable]
                self.gainfields += ['']
            self.logger.info('Done.')
        #
        # Correct for gaincurve and antenna efficiencies
        #
        if self.gaincurve:
            caltable = 'gaincurve.cal'
            self.logger.info('Calculating gain curve and antenna '
                             'efficiencies...')
            if os.path.isdir(caltable):
                casa.rmtables(caltable)
            casa.gencal(vis=self.vis, caltable=caltable,
                        caltype='gceff')
            if os.path.isdir(caltable):
                self.gaintables += [caltable]
                self.gainfields += ['']
            self.logger.info('Done.')
        #
        # Correct for atmospheric opacity
        #
        if self.opacity:
            caltable = 'opacity.cal'
            self.logger.info('Calculating opacity correction...')
            if os.path.isdir(caltable):
                casa.rmtables(caltable)
            myTau = casa.plotweather(vis=self.vis, doPlot=False)
            allspws = ','.join([str(spw) for spw in range(len(myTau))])
            casa.gencal(vis=self.vis, caltable=caltable, caltype='opac',
                        parameter=myTau, spw=allspws)
            if os.path.isdir(caltable):
                self.gaintables += [caltable]
                self.gainfields += ['']
            self.logger.info('Done.')
        #
        # pre-bandpass calibration delay calibration on primary
        # calibrators (linear slope in phase vs frequency)
        #
        caltable = 'delays.Kcal'
        self.logger.info('Calculating delay calibration table for '
                         'primary calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.pri_cals)
        gaintables = self.gaintables
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     refant=self.refant, gaintype='K', minblperant=1,
                     parang=self.parang, gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with delay calibration')
            raise ValueError('Problem with delay calibration!')
        self.delays = caltable
        self.logger.info('Done.')
        #
        # pre-bandpass integration timescale phase calibration
        # (phase vs time)
        #
        caltable = 'phase_int.Gcal0'
        self.logger.info('Calculating phase calibration table on '
                         'integration timescales for primary '
                         'calibrators...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.pri_cals)
        gaintables = self.gaintables + [self.delays]
        casa.gaincal(vis=self.vis, caltable=caltable, field=field,
                     solint='int', calmode='p', refant=self.refant,
                     gaintype='G', minsnr=2.0, minblperant=1,
                     parang=self.parang, gaintable=gaintables)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with integration-timescale '
                                 'phase calibration')
            raise ValueError('Problem with integration-timescale '
                             'phase calibration!')
        self.phase_int = caltable
        self.logger.info('Done.')
        #
        # Remaining calibration tables
        #
        self.calibration_tables()
        #
        # Polarization calibration
        #
        if self.calpol:
            #
            # Re-do bandpass and gain calibrations with polarization
            # leakage tables
            #
            self.calibration_tables_post_polcal()
        #
        # set the flux scale
        #
        caltable = 'flux.cal'
        self.logger.info('Calculating the flux calibration table...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        field = ','.join(self.flux_cals)
        casa.fluxscale(vis=self.vis, caltable=self.apcal_scan,
                       fluxtable=caltable, reference=field,
                       incremental=True)
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with flux calibration')
            raise ValueError('Problem with flux calibration!')
        self.flux = caltable
        self.logger.info('Done.')
        #
        # Apply calibration solutions to calibrators
        #
        self.logger.info('Applying calibration tables to all '
                         'calibrators...')
        if self.calpol:
            gaintables = self.gaintables + [
                self.delays, self.bandpass, self.phase_int,
                self.apcal_scan, self.polcal_scan, self.flux]
        else:
            gaintables = self.gaintables + [
                self.delays, self.bandpass, self.phase_int,
                self.apcal_scan, self.flux]
        for field in self.pri_cals+self.sec_cals:
            if self.calpol:
                gainfields = self.gainfields + [
                    '', '', field, field, '', field]
            else:
                gainfields = self.gainfields + [
                    '', '', field, field, field]
            casa.applycal(vis=self.vis, field=field, calwt=self.calwt,
                          parang=self.parang, gaintable=gaintables,
                          gainfield=gainfields, flagbackup=False)
        self.logger.info('Done.')
        self.save_flags('calibrate')
        #
        # Generate calibration plots
        #
        field = ','.join(self.pri_cals)
        casa.plotms(vis=self.bandpass, xaxis='channel',
                    yaxis='amplitude', field=field, iteraxis='spw',
                    coloraxis='antenna1',
                    title=self.bandpass.replace('_', '\_'),
                    plotfile='plotcal_figures/0_bandpass.png',
                    overwrite=True, showgui=False, exprange='all')
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.plotms(vis=self.phase_int, xaxis='time', yaxis='phase',
                    field=field, iteraxis='spw',
                    coloraxis='antenna1',
                    title=self.phase_int.replace('_', '\_'),
                    plotfile='plotcal_figures/1_phase_int.png',
                    overwrite=True, showgui=False, exprange='all')
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.plotms(vis=self.phase_scan, xaxis='time', yaxis='phase',
                    field=field, iteraxis='spw',
                    coloraxis='antenna1',
                    title=self.phase_scan.replace('_', '\_'),
                    plotfile='plotcal_figures/2_phase_scan.png',
                    overwrite=True, showgui=False, exprange='all')
        field = ','.join(self.pri_cals+self.sec_cals)
        casa.plotms(vis=self.apcal_scan, xaxis='time',
                    yaxis='amplitude', field=field, iteraxis='spw',
                    coloraxis='antenna1',
                    title=self.apcal_scan.replace('_', '\_'),
                    plotfile='plotcal_figures/3_apcal_scan.png',
                    overwrite=True, showgui=False, exprange='all')
        #
        # Generate PDF of plotcal figures
        #
        self.logger.info('Generating tex document...')
        iplot = 0
        with open('plotcal_figures.tex','w') as fout:
            fout.write(r'\documentclass{article}'+'\n')
            fout.write(r'\usepackage{graphicx,subfig}'+'\n')
            fout.write(r'\usepackage[margin=0.1cm]{geometry}'+'\n')
            fout.write(r'\begin{document}'+'\n')
            fout.write(r'\begin{figure}'+'\n')
            fout.write(r'\centering'+'\n')
            fnames = glob.glob('plotcal_figures/*.png')
            fnames = natural_sort(fnames)
            for fname in fnames:
                if iplot > 0 and iplot % 6 == 0:
                    fout.write(r'\end{figure}'+'\n')
                    fout.write(r'\clearpage'+'\n')
                    fout.write(r'\begin{figure}'+'\n')
                    fout.write(r'\centering'+'\n')
                elif iplot > 0 and iplot % 2 == 0:
                    fout.write(r'\end{figure}'+'\n')
                    fout.write(r'\begin{figure}'+'\n')
                    fout.write(r'\centering'+'\n')
                fout.write(r'\includegraphics[width=0.45\textwidth]{'+fname+'}'+'\n')
                iplot += 1
            fout.write(r'\end{figure}'+'\n')
            fout.write(r'\end{document}'+'\n')
        os.system('pdflatex -interaction=batchmode plotcal_figures.tex')
        self.logger.info('Done.')

    def calibrate_sci_targets(self):
        """
        Apply calibration solutions to science targets.

        Inputs: Nothing

        Returns: Nothing
        """
        if self.calpol:
            gaintables = self.gaintables + [
                self.delays, self.bandpass, self.phase_scan,
                self.apcal_scan, self.polcal_scan, self.flux]
            gainfields = self.gainfields + [
                '', '', 'nearest', 'nearest', '', 'nearest']
        else:
            gaintables = self.gaintables + [
                self.delays, self.bandpass, self.phase_scan,
                self.apcal_scan, self.flux]
            gainfields = self.gainfields + [
                '', '', 'nearest', 'nearest', 'nearest']
        for field in self.sci_targets:
            self.logger.info('Applying calibration solutions to %s',
                             field)
            casa.applycal(vis=self.vis, field=field, calwt=self.calwt,
                          gaintable=gaintables, gainfield=gainfields,
                          parang=self.parang, flagbackup=False)
        self.save_flags('calibrate')

    def auto_flag_sci_targets(self):
        """
        Perform automatic flagging of calibrated science targets
        using rflag.

        Inputs: Nothing

        Returns: Nothing
        """
        field = ','.join(self.sci_targets)
        datacolumn = 'corrected'
        self.logger.info('Running rflag on all correlations...')
        casa.flagdata(vis=self.vis, mode='rflag', field=field,
                      flagbackup=False, datacolumn=datacolumn,
                      extendflags=False)
        self.save_flags('autoflag')

    def sci_target_plots(self):
        """
        Generate science target diagnostic plots.

        Inputs: Nothing

        Returns: Nothing
        """
        self.logger.info('Generating plots for manual inspection...')
        plotnum = 0
        plots = []
        corr = self.config.get('Polarization', 'Polarization')
        datacolumn = 'corrected'
        for field in self.sci_targets:
            #
            # Amplitude vs UV-distance (in wavelength units)
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'science_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='uvwave', yaxis='amp',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', coloraxis='baseline',
                        correlation=corr, title=title,
                        plotfile=plotfile, overwrite=True,
                        showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'uvwave', 'yaxis':'amp',
                          'avgtime':'', 'avgchannel':''})
            plotnum += 1
            #
            # Amplitude vs Time
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'science_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='time', yaxis='amp',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', avgchannel='1e7',
                        coloraxis='baseline', correlation=corr,
                        overwrite=True, showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'time', 'yaxis':'amp',
                          'avgtime':'', 'avgchannel':'1e7'})
            plotnum += 1
            #
            # Amplitude vs Channel
            #
            title = 'PlotID: {0} Field: {1}'.format(plotnum, field)
            plotfile = 'science_figures/{0}.png'.format(plotnum)
            casa.plotms(vis=self.vis, xaxis='channel', yaxis='amp',
                        field=field, ydatacolumn=datacolumn,
                        iteraxis='spw', avgtime='1e7',
                        coloraxis='baseline', correlation=corr,
                        showgui=False, exprange='all')
            plots.append({'field':field,
                          'xaxis':'channel', 'yaxis':'amp',
                          'avgtime':'1e7', 'avgchannel':''})
            plotnum += 1
        self.logger.info('Done.')
        #
        # Generate PDF to display plots
        #
        self.logger.info('Generating tex document...')
        num_plots = plotnum
        iplot = 0
        with open('science_figures.tex', 'w') as fout:
            fout.write(r'\documentclass{article}'+'\n')
            fout.write(r'\usepackage{graphicx}'+'\n')
            fout.write(r'\usepackage[margin=0.1cm]{geometry}'+'\n')
            fout.write(r'\begin{document}'+'\n')
            fout.write(r'\begin{figure}'+'\n')
            fout.write(r'\centering'+'\n')
            for plotnum in range(num_plots):
                fnames = glob.glob('science_figures/{0}_*.png'.
                                   format(plotnum))
                fnames = natural_sort(fnames)
                for fname in fnames:
                    if iplot > 0 and iplot % 6 == 0:
                        fout.write(r'\end{figure}'+'\n')
                        fout.write(r'\clearpage'+'\n')
                        fout.write(r'\begin{figure}'+'\n')
                        fout.write(r'\centering'+'\n')
                    elif iplot > 0 and iplot % 2 == 0:
                        fout.write(r'\end{figure}'+'\n')
                        fout.write(r'\begin{figure}'+'\n')
                        fout.write(r'\centering'+'\n')
                    fout.write(r'\includegraphics[width=0.45\textwidth]{'+fname+'}\n')
                    iplot += 1
            fout.write(r'\end{figure}'+'\n')
            fout.write(r'\end{document}'+'\n')
        os.system('pdflatex -interaction=batchmode science_figures.tex')
        self.logger.info('Done.')
        #
        # Save the plot list to the pickle object
        #
        self.logger.info('Saving plot list to pickle...')
        with open('science_plots.pkl', 'w') as fout:
            pickle.dump(plots, fout)
        self.logger.info('Done.')

    def manual_flag_sci_targets(self):
        """
        Interactively plot and flag the science targets.

        Inputs: Nothing

        Returns: Nothing
        """
        #
        # Read plot list from pickle object
        #
        self.logger.info('Reading plot list from pickle...')
        with open('science_plots.pkl', 'r') as fin:
            plots = pickle.load(fin)
        num_plots = len(plots)
        self.logger.info('Done.')
        #
        # Prompt user with menu
        #
        corr = self.config.get('Polarization', 'Polarization')
        self.logger.info('Please inspect science_plots.pdf then perform manual flagging.')
        while True:
            print('f - flag some data')
            print('plot id number - generate interactive version of plot with this id')
            print('quit - end this flagging session')
            answer = input()
            #
            # Flag some data
            #
            if answer.lower() == 'f':
                self.flag(self.sci_targets)
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
                    print('Invalid Plot ID')
                    continue
                if plotid >= num_plots:
                    print('Invalid Plot ID')
                    continue
                title = 'PlotID: {0} Field: {1}'.format(
                    plotid, plots[plotid]['field'])
                casa.plotms(vis=self.vis,
                            xaxis=plots[plotid]['xaxis'],
                            yaxis=plots[plotid]['yaxis'],
                            field=plots[plotid]['field'],
                            ydatacolumn='corrected', iteraxis='spw',
                            title=title, coloraxis='baseline',
                            correlation=corr,
                            avgchannel=plots[plotid]['avgchannel'],
                            avgtime=plots[plotid]['avgtime'])
        #
        # Save the flags
        #
        self.save_flags('manualflag')

    def split_fields(self):
        """
        Split calibrated fields into measurement sets with naming:
        {field_name}_calibrated.ms

        Inputs: Nothing

        Returns: Nothing
        """
        for field in self.all_fields:
            outputvis = '{0}_calibrated.ms'.format(field)
            self.logger.info('Splitting %s to %s', field, outputvis)
            casa.split(vis=self.vis, outputvis=outputvis, field=field,
                       keepflags=False)
        self.logger.info('Done!')

def main(vis, config_file, shadow_tolerance=0.0,
         quack_interval=10.0, antpos=True, gaincurve=True,
         opacity=True, calpol=False, calwt=True, auto=''):
    """
    Run the calibration pipeline

    Inputs:
      vis :: string
        The masurement set
      config_file :: string
        The filename of the configuration file for this project
      shadow_tolerance :: scalar
        The overlap tolerance used for shadow flagging. Flag
        any data with projected baseline separation less than
        r_1 + r_2 - shadow_tolerance
        where r_1 and r_2 are the radii of the antennas.
      quack_interval :: scalar
        The amount of time in seconds to flag at the beginning
        of each scan.
      antpos :: boolean
        if True, compute antenna position corrections (only for VLA)
      gaincurve :: boolean
        if True, compute gain curve and antenna efficiency corrections
        (only for VLA)
      opacity :: boolean
        if True, compute opacity corrections
      calpol :: boolean
        if True, calibrate polarization
      calwt :: boolean
        if True, apply calibration weights to data
      auto :: string
        comma separated string of menu options to automatically
        perform

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger('main')
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
    logger.info('Reading configuration file %s', config_file)
    config.read(config_file)
    logger.info('Done.')
    #
    # Initialize Calibration object
    #
    calib = Calibration(
        vis, logger, config, shadow_tolerance=shadow_tolerance,
        quack_interval=quack_interval, antpos=antpos,
        gaincurve=gaincurve, opacity=opacity, calpol=calpol,
        calwt=calwt)
    #
    # Prompt the user with a menu, or automatically execute
    #
    auto_items = auto.split(',')
    auto_ind = 0
    while True:
        if not auto:
            print('0. Preliminary flags (config file, etc.)')
            print('1. Auto-flag calibrator fields')
            print('2. Generate plotms figures for calibrator fields')
            print('3. Manually flag calibrator fields')
            print('4. Calculate and apply calibration solutions to calibrator fields')
            print('5. Apply calibration solutions to science fields')
            print('6. Auto-flag science fields')
            print('7. Generate plotms figures for science fields')
            print('8. Manually flag science fields')
            print('9. Split calibrated fields')
            print('q [quit]')
            answer = input('> ')
        else:
            answer = auto_items[auto_ind]
            auto_ind += 1
        if answer == '0':
            calib.preliminary_flagging()
        elif answer == '1':
            calib.auto_flag_calibrators()
        elif answer == '2':
            calib.calibrator_plots()
        elif answer == '3':
            calib.manual_flag_calibrators()
        elif answer == '4':
            calib.calibrate_calibrators()
        elif answer == '5':
            calib.calibrate_sci_targets()
        elif answer == '6':
            calib.auto_flag_sci_targets()
        elif answer == '7':
            calib.sci_target_plots()
        elif answer == '8':
            calib.manual_flag_sci_targets()
        elif answer == '9':
            calib.split_fields()
        elif answer.lower() == 'q' or answer.lower() == 'quit':
            break
        else:
            print('Input not recognized.')
        if auto_ind >= len(auto_items):
            break
