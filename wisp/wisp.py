"""
wisp.py - Wenger Interferometry Software Package

Copyright(C) 2018-2021 by
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
Trey V. Wenger - April 2021 - v3.0
    Improve code readability.
"""

import os
import logging
import logging.config
import ConfigParser

from .calibration import Calibration, apply_calibration, generate_tables
from .flagging import auto_flag, manual_flag, preliminary_flagging
from .calibration_plots import visibility_plots, plotcal_plots

from .imaging import Imaging
from .mfs_images import (
    mfs_dirty_cont,
    mfs_clean_cont,
    mfs_dirty_spws,
    mfs_clean_spws,
)
from .channel_images import channel_dirty_spws, channel_clean_spws
from .imaging_plots import contplot, lineplot

import __main__ as casa

# Catch raw_input in python 3
try:
    input = raw_input
except NameError:
    raw_input = input

# load logging configuration file
logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "logging.conf"))

__version__ = "3.1"


def calibrate(
    vis,
    config_file,
    refant,
    shadow_tolerance=0.0,
    quack_interval=0.0,
    antpos=True,
    gaincurve=True,
    opacity=True,
    calpol=False,
    calwt=True,
    solint="int",
    auto="",
):
    """
    Run the calibration pipeline

    Inputs:
        vis :: string
            The masurement set
        config_file :: string
            The filename of the configuration file for this project
        refant :: string
            Reference antenna
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
        solint :: string
            The solution interval for short timescale phase corrections.
            Default is 'int' for integration. (e.g., solint='10s')
        auto :: string
            comma separated string of menu options to automatically
            perform

    Returns: Nothing
    """

    # start logger
    logger = logging.getLogger("main")

    # Check inputs
    if not os.path.isdir(vis):
        logger.critical("Measurement set not found!")
        raise ValueError("Measurement set not found!")
    if not os.path.exists(config_file):
        logger.critical("Configuration file not found")
        raise ValueError("Configuration file not found!")

    # load configuration file
    config = ConfigParser.ConfigParser()
    logger.info("Reading configuration file %s", config_file)
    config.read(config_file)
    logger.info("Done.")

    # Initialize Calibration object
    cal = Calibration(
        casa,
        vis,
        logger,
        config,
        refant,
        shadow_tolerance=shadow_tolerance,
        quack_interval=quack_interval,
        antpos=antpos,
        gaincurve=gaincurve,
        opacity=opacity,
        calpol=calpol,
        calwt=calwt,
        solint=solint,
    )

    # Prompt the user with a menu, or automatically execute
    auto_items = auto.split(",")
    auto_ind = 0
    while True:
        if not auto:
            print("0. Preliminary flags (config file, etc.)")
            print("1. Auto-flag calibrator fields")
            print("2. Generate plotms figures for calibrator fields")
            print("3. Manually flag calibrator fields")
            print(
                "4. Calculate and apply calibration solutions to calibrator " "fields"
            )
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
        if answer == "0":
            preliminary_flagging(cal)
        elif answer == "1":
            auto_flag(cal, cal.calibrators, extend=True)
        elif answer == "2":
            visibility_plots(cal, "calibrator")
        elif answer == "3":
            manual_flag(cal, "calibrator")
        elif answer == "4":
            generate_tables(cal)
            plotcal_plots(cal)
            apply_calibration(cal, "calibrator")
        elif answer == "5":
            apply_calibration(cal, "science")
        elif answer == "6":
            auto_flag(cal, cal.sci_targets, extend=False)
        elif answer == "7":
            visibility_plots(cal, "science")
        elif answer == "8":
            manual_flag(cal, "science")
        elif answer == "9":
            cal.split_fields()
        elif answer.lower() == "q" or answer.lower() == "quit":
            break
        else:
            print("Input not recognized.")
        if auto_ind >= len(auto_items):
            break


def imaging(
    vis,
    field,
    config_file,
    outdir=".",
    stokes="I",
    spws="",
    uvrange="",
    uvtaper=False,
    outertaper="",
    imsize=None,
    phasecenter=None,
    interactive=False,
    savemodel=None,
    parallel=False,
    auto="",
):
    """
    Generate and clean images

    Inputs:
      vis :: string
        The measurement set containing all data for field
      field :: string
        The field name to image
      config_file :: string
        filename of the configuration file for this project
      outdir :: string
        The directory where to save the results. This is useful if
        you plan to make several sets of images (i.e. self-calibration)
      stokes :: string
        The Stokes parameters we're imaging. e.g. 'I' or 'IQUV'
      spws :: string
        comma-separated list of spws and channels to clean
        if empty, clean all spws, all channels
      uvrange :: string
        Selection on UV-range
      uvtaper :: boolean
        if True, apply UV tapering
      outertaper :: string
        Tapering FWHM
      imsize :: list of integers
        Image size. If None, use config file.
      phasecenter :: string
        If not None, use this phase center
      interactive :: boolean
        if True, interactively clean
      savemodel :: string
        if not none, save individual MFS images of each spectral
        window to the model column of the measurement set for
        self-calibration. This can only be done with stokes='I'.
        if savemodel == 'light': save the model after lightniter
        if savemodel == 'clean': save the model after niter
      parallel :: boolean
        if True, run parallel TCLEAN
        N.B. CASA must be started in MPI mode (via mpicasa)
      auto :: string
        if not an empty string, it is a comma separated
        list of menu items to perform, i.e. auto='0,1,4,5,6'

    Returns: Nothing
    """
    #
    # start logger
    #
    logger = logging.getLogger("main")
    #
    # Check inputs
    #
    if not os.path.isdir(vis):
        logger.critical("Measurement set not found!")
        raise ValueError("Measurement set not found!")
    if not os.path.exists(config_file):
        logger.critical("Configuration file not found")
        raise ValueError("Configuration file not found!")
    if savemodel is not None and stokes != "I":
        logger.critical("Can only save visibility model with Stokes I")
        raise ValueError("Can only save visibility model with Stokes I")
    #
    # load configuration file
    #
    config = ConfigParser.ConfigParser()
    logger.info("Reading configuration file {0}".format(config_file))
    config.read(config_file)
    logger.info("Done.")
    #
    # Initialize Imaging object
    #
    img = Imaging(
        vis,
        field,
        logger,
        config,
        outdir=outdir,
        uvtaper=uvtaper,
        outertaper=outertaper,
        spws=spws,
        uvrange=uvrange,
        stokes=stokes,
        imsize=imsize,
        phasecenter=phasecenter,
        savemodel=savemodel,
        interactive=interactive,
        parallel=parallel,
    )
    #
    # Prompt the user with a menu for each option, or auto-do them
    #
    auto_items = auto.split(",")
    auto_ind = 0
    while True:
        if not auto:
            print(
                "0. Dirty image combined continuum spws "
                "(MFS; multi-term; multi-scale)"
            )
            print("1. Clean combined continuum spws " "(MFS; multi-term; multi-scale)")
            print("2. Dirty image each continuum spw (MFS; multi-scale)")
            print("3. Clean each continuum spw (MFS; multi-scale)")
            print("4. Dirty image each continuum spw (channel; multi-scale)")
            print("5. Clean each continuum spw (channel; multi-scale)")
            print("6. Dirty image each line spw (MFS; multi-scale)")
            print("7. Clean each line spw (MFS; multi-scale)")
            print("8. Dirty image each line spw (channel; multi-scale)")
            print("9. Clean each line spw (channel; multi-scale)")
            print("10. Generate continuum diagnostic plots")
            print("11. Generate spectral line diagnostic plots")
            print("q [quit]")
            answer = raw_input("> ")
        else:
            answer = auto_items[auto_ind]
            auto_ind += 1
        if answer == "0":
            mfs_dirty_cont(img)
        elif answer == "1":
            mfs_clean_cont(img)
        elif answer == "2":
            mfs_dirty_spws(img, img.cont_spws, img.cont_chans)
        elif answer == "3":
            mfs_clean_spws(img, img.cont_spws, img.cont_chans, "cont")
        elif answer == "4":
            channel_dirty_spws(img, img.cont_spws, "cont")
        elif answer == "5":
            channel_clean_spws(img, img.cont_spws, "cont")
        elif answer == "6":
            mfs_dirty_spws(img, img.line_spws, img.line_chans)
        elif answer == "7":
            mfs_clean_spws(img, img.line_spws, img.line_chans, "line")
        elif answer == "8":
            channel_dirty_spws(img, img.line_spws, "line")
        elif answer == "9":
            channel_clean_spws(img, img.line_spws, "line")
        elif answer == "10":
            contplot(img)
        elif answer == "11":
            lineplot(img)
        elif answer.lower() == "q" or answer.lower() == "quit":
            break
        else:
            print("Input not recognized.")
        if auto_ind >= len(auto_items):
            break
