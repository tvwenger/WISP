"""
reweight.py - WISP re-compute data weights

Re-compute data weights using CASA task STATWT to estimate noise
from visibilities.

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
"""

import __main__ as casa
import os
import time
import numpy as np
import logging
import logging.config
import ConfigParser

__version__ = "1.0"

# load logging configuration file
logging.config.fileConfig('logging.conf')

def setup(vis='',config=None):
    """
    Perform setup tasks: find line and continuum spectral windows
                         get clean parameters

    Inputs:
      vis :: string
        The measurement set
      config :: a ConfigParser object
        The ConfigParser object for this project

    Returns: my_cont_spws, my_line_spws, chan_offset, chan_width
      my_cont_spws :: string
        comma-separated string of continuum spws
      my_line_spws :: string
        comma-separated string of line spws
      chan_offset :: 1-D array of integers
        The channel offsets relative to last line spectral window
      chan_width :: integer
        The total number of channels to un-flag
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
    # Get continuum and line spws from configuration file
    #
    my_cont_spws = config.get("Spectral Windows","Continuum")
    my_line_spws = config.get("Spectral Windows","Line")
    logger.info("Found continuum spws: {0}".format(my_cont_spws))
    logger.info("Found line spws: {0}".format(my_line_spws))
    #
    # Get Unflag parameters from configuration file
    #
    chan_offset = np.array([int(c) for c in config.get("Unflag","offset").split(',')])
    chan_width = config.getint("Unflag","width")
    return (my_cont_spws,my_line_spws,chan_offset,chan_width)

def main(field,vis='',config_file=''):
    """
    Re-compute data weights using STATWT. Use "unflag" parameters
    from configuration file to exclude spectral lines.

    Inputs:
      field :: string
        The field name to clean
      vis :: string
        The measurement set containing all data for field
      config_file :: string
        filename of the configuration file for this project

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Check inputs
    #
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
    my_cont_spws,my_line_spws,chan_offset,chan_width = \
      setup(vis=vis,config=config)
    #
    # Get spw we are looking at
    #
    print("What spectral window are you using to measure the line channel?")
    print("<CR> to indicate no visible spectral line")
    spw = raw_input("> ")
    if len(spw.strip()) > 0:
        spw_ind = my_line_spws.split(',').index(spw)
        #
        # Get channel of spectral line in highest-freq spectral window
        #
        print("What is the channel of the spectral line center in that spectral wnidow?")
        chan = int(raw_input("> "))
        #
        # Chan offset is relative to last spw, so we need to adjust
        # the offsets for this spw
        #
        chan_offset = chan_offset + chan_offset[-1] - chan_offset[spw_ind]
        #
        # Backup flags
        #
        logger.info("Saving flag state...")
        casa.flagmanager(vis=vis,mode='save',versionname='flags_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
        logger.info("Done.")
        #
        # Run statwt on continuum spws
        #
        logger.info("Running statwt on continuum spws")
        casa.statwt(vis=vis,spw=my_cont_spws,datacolumn='data',flagbackup=False)
        logger.info("Done.")
        #
        # Run statwt on line spws
        #
        for spw, offset in zip(my_line_spws.split(','),chan_offset):
            start = chan + int(offset) - int(chan_width/2)
            end = chan + int(offset) + int(chan_width/2)
            logger.info("Running statwt on spw {0} excluding channels {1} to {2}".format(spw,start,end))
            casa.statwt(vis=vis,spw=spw,datacolumn='data',flagbackup=False,
                        excludechans='{0}:{1}~{2}'.format(spw,start,end))
            logger.info("Done.")
    else:
        #
        # Backup flags
        #
        logger.info("Saving flag state...")
        casa.flagmanager(vis=vis,mode='save',versionname='flags_{0}'.format(time.strftime('%Y%m%d%H%M%S',time.gmtime())))
        logger.info("Done.")
        #
        # Execute statwt on all spws
        #
        logger.info("Running statwt on all spws")
        casa.statwt(vis=vis,datacolumn='data',flagbackup=False)
        logger.info("Done.")
