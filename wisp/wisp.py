"""
wisp.py - Wenger Interferometry Software Package

Copyright(C) 2018-2020 by
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
Trey V. Wenger August 2019 - V2.1
    Improve code readability.
"""

import os
import logging
import logging.config
import ConfigParser

from .calibration import Calibration, apply_calibration, generate_tables
from .flagging import auto_flag, manual_flag, preliminary_flagging
from .plots import visibility_plots

import __main__ as casa

# load logging configuration file
logging.config.fileConfig(
    os.path.join(os.path.dirname(__file__), "logging.conf")
)


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
                "4. Calculate and apply calibration solutions to calibrator "
                "fields"
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
