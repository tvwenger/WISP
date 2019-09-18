"""
selfcal.py - WISP Self-calibration Pipeline

Self-calibrate a measurement set.

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
Trey V. Wenger August 2019 - V2.0
   Initial version in v2.0 of WISP.
"""

import os
import glob
import re
import time

import ConfigParser
import logging
import logging.config

import __main__ as casa

__version__ = '2.0'

# load logging configuration file
logging.config.fileConfig('logging.conf')

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

class SelfCalibration:
    """
    The SelfCalibration object handles the self-calibration steps for
    a measurement set.
    """

    def __init__(self, vis, refant, logger, config):
        """
        Create a new SelfCalibration object. Create a listobs file
        and generate plotcal directory.

        Inputs:
          vis :: string
            The masurement set
          refant :: string
            The reference antenna to use
          logger :: logging.Logger object
            The logging object we're using
          config :: config.ConfigParser object
            The config parser we're using

        Returns: self_calibration
          self_calibration :: selfcal.SelfCalibration object
            a new SelfCalibration object
        """
        self.vis = vis
        self.refant = refant
        self.logger = logger
        self.config = config
        self.logger.info('Initializing SelfCalibration object.')
        #
        # Generate listobs file
        #
        listfile = 'listobs.txt'
        if not os.path.isfile(listfile):
            self.logger.info('Generating listobs file...')
            casa.listobs(vis=self.vis, listfile=listfile)
            self.logger.info('Done.')
        #
        # get field names
        #
        self.logger.info('Looking for field names...')
        all_fields = casa.vishead(
            vis=self.vis, mode='get', hdkey='field')[0]
        all_fields = list(set(all_fields))
        self.logger.info('Found fields: %s',
                         ', '.join(all_fields))
        if len(all_fields) > 1:
            self.logger.critical('There can only be one field in MS!')
            raise ValueError('There can only be one field in MS!')
        self.field = all_fields[0]
        #
        # Determine if we are doing parallactic angle correction
        #
        self.parang = self.config.getboolean(
            'Polarization', 'Parallactic Angle Correction')
        #
        # create directories for figures
        #
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

    def selfcal_phase(self):
        """
        Generate phase calibration table.

        Inputs: Nothing

        Returns: Nothing
        """
        caltable = 'phase_scan.Gcal'
        self.logger.info('Calculating the phase calibration table on '
                         'scan timescales...')
        if os.path.isdir(caltable):
            casa.rmtables(caltable)
        casa.gaincal(vis=self.vis, caltable=caltable, field=self.field,
                     solint='inf', calmode='p', refant=self.refant,
                     gaintype='G', minsnr=2.0, minblperant=1,
                     parang=self.parang, gaintable=[])
        if not os.path.isdir(caltable):
            self.logger.critical('Problem with scan-timescale phase '
                                 'calibration')
            raise ValueError('Problem with scan-timescale phase '
                             'calibration!')
        self.logger.info('Done.')

    def calibrate(self):
        """
        Calculate calibration solutions and apply solutions to
        the self-calibrators.

        Inputs: Nothing

        Returns: Nothing
        """
        self.logger.info('Applying calibration tables...')
        for field in [self.field]:
            # N.B. No calwt for self-calibration
            casa.applycal(vis=self.vis, field=field, calwt=False,
                          parang=self.parang,
                          gaintable=['phase_scan.Gcal'],
                          gainfield=[self.field],
                          flagbackup=False)
        self.logger.info('Done.')
        self.save_flags('post self-calibrate')
        #
        # Remove previous figures
        #
        fnames = glob.glob('plotcal_figures/*.png')
        for fname in fnames:
            os.remove(fname)
        #
        # Generate calibration plots
        #
        casa.plotms(vis='phase_scan.Gcal', xaxis='time', yaxis='phase',
                    field=self.field, iteraxis='spw',
                    coloraxis='antenna1',
                    title='phase_scan.Gcal'.replace('_', '\_'),
                    plotfile='plotcal_figures/0_phase_scan.png',
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

def main(vis, refant, config_file):
    """
    Run the self-calibration pipeline

    Inputs:
      vis :: string
        The masurement set
      refant :: string
        The reference antenna to use
      config_file :: string
        The filename of the configuration file for this project

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
    # Initialize SelfCalibration object
    #
    calib = SelfCalibration(
        vis, refant, logger, config)
    #
    # Save initial flags
    #
    calib.save_flags('pre self-calibrate')
    #
    # Generate phase self-calibration table and apply
    #
    calib.selfcal_phase()
    calib.calibrate()
