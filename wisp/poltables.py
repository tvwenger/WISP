"""
poltables.py - Generate polarization calibration tables

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


def crosshand_delays_table(cal, casa):
    """
    Derive cross-hand delay calibration table using a bright,
    polarized calibrator.

    Inputs:
        cal :: Calibration object
            The calibration object
        casa :: CASA namespace
            CASA namespace

    Returns: Nothing
    """
    field = ",".join(cal.pol_leak_cals)
    initial_tables = cal.initial_tables()

    caltable = "crosshand_delays.Kcal"
    gaintables = initial_tables + [
        cal.table["delays"],
        cal.table["phase_int"],
        cal.table["apcal_scan"],
    ]
    cal.logger.info(
        "Calculating cross-hand delay calibration table for polarized "
        "leakage calibrators..."
    )
    if os.path.isdir(caltable):
        casa.rmtables(caltable)
    casa.gaincal(
        vis=cal.vis,
        caltable=caltable,
        field=field,
        refant=cal.refant,
        solint="inf",
        combine="scan,spw",
        minblperant=1,
        gaintype="KCROSS",
        parang=cal.calpol,
        gaintable=gaintables,
    )
    if not os.path.isdir(caltable):
        cal.logger.critical("Problem with cross-hand delay calibration")
        raise ValueError("Problem with cross-hand delay calibration!")
    cal.tables["crosshand_delays_scan_spw"] = caltable
    cal.logger.info("Done.")


def polleak_table(cal, casa):
    """
    Derive polarization leakage calibration table.

    Inputs:
        cal :: Calibration object
            The calibration object
        casa :: CASA namespace
            CASA namespace

    Returns: Nothing
    """
    field = ",".join(cal.pol_leak_cals)
    initial_tables = cal.initial_tables()

    caltable = "polleak_scan.Dcal0"
    gaintables = initial_tables + [
        cal.tables["delays"],
        cal.tables["bandpass"],
        cal.tables["phase_int"],
        cal.tables["apcal_scan"],
    ]
    if "crosshand_delays_scan_spw" in cal.tables:
        gaintables = gaintables + [cal.tables["crosshand_delays_scan_spw"]]
    cal.logger.info(
        "Calculating the polarization leakage calibration table on scan "
        "timescales for polarization leakage calibrators ..."
    )
    if os.path.isdir(caltable):
        casa.rmtables(caltable)
    casa.polcal(
        vis=cal.vis,
        caltable=caltable,
        field=field,
        solint="inf",
        combine="scan",
        poltype="Df",
        refant=cal.refant,
        minsnr=2.0,
        minblperant=1,
        gaintable=gaintables,
    )
    if not os.path.isdir(caltable):
        cal.logger.critical("Problem with polarization leakage calibration")
        raise ValueError("Problem with polarization leakage calibration")
    cal.tables["polleak_scan"] = caltable
    cal.logger.info("Done.")


def polangle_table(cal, casa):
    """
    Derive polarization angle calibration table.

    Inputs:
        cal :: Calibration object
            The calibration object
        casa :: CASA namespace
            CASA namespace

    Returns: Nothing
    """
    field = ",".join(cal.pol_angle_cals)
    initial_tables = cal.initial_tables()

    caltable = "polangle_scan.Xcal"
    gaintables = initial_tables + [
        cal.tables["delays"],
        cal.tables["bandpass"],
        cal.tables["phase_int"],
        cal.tables["apcal_scan"],
    ]
    if "crosshand_delays_scan_spw" in cal.tables:
        gaintables = gaintables + [cal.tables["crosshand_delays_scan_spw"]]
    gaintables = gaintables + [cal.tables["polleak_scan"]]
    cal.logger.info(
        "Calculating the polarization angle calibration table on scan "
        "timescales for polarization angle calibrators ..."
    )
    if os.path.isdir(caltable):
        casa.rmtables(caltable)
    casa.polcal(
        vis=cal.vis,
        caltable=caltable,
        field=field,
        solint="inf",
        combine="scan",
        poltype="Xf",
        refant=cal.refant,
        minsnr=2.0,
        minblperant=1,
        gaintable=gaintables,
    )
    if not os.path.isdir(caltable):
        cal.logger.critical("Problem with polarization angle calibration")
        raise ValueError("Problem with polarization angle calibration")
    cal.tables["polangle_scan"] = caltable
    cal.logger.info("Done.")
