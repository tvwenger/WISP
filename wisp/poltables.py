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


def crosshand_delays_table(cal):
    """
    Derive cross-hand delay calibration table using a bright,
    polarized calibrator.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: Nothing
    """
    fields = cal.pol_leak_cals

    cal.logger.info(
        "Calculating cross-hand delay calibration table for polarized "
        "leakage calibrators..."
    )
    if os.path.isdir(cal.tables["crosshand_delays"]):
        cal.casa.rmtables(cal.tables["crosshand_delays"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables(
            "crosshand_delays", field
        )
        append = os.path.exists(cal.tables["crosshand_delays"])
        cal.casa.gaincal(
            vis=cal.vis,
            caltable=cal.tables["crosshand_delays"],
            field=field,
            refant=cal.refant,
            solint="inf",
            combine="scan,spw",
            minblperant=1,
            gaintype="KCROSS",
            parang=cal.calpol,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
        )
    if not os.path.isdir(cal.tables["crosshand_delays"]):
        cal.logger.critical("Problem with cross-hand delay calibration")
        raise ValueError("Problem with cross-hand delay calibration!")
    cal.logger.info("Done.")


def polleak_table(cal):
    """
    Derive polarization leakage calibration table.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: Nothing
    """
    fields = cal.pol_leak_cals

    cal.logger.info(
        "Calculating the polarization leakage calibration table on scan "
        "timescales for polarization leakage calibrators ..."
    )
    if os.path.isdir(cal.tables["polleak"]):
        cal.casa.rmtables(cal.tables["polleak"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("polleak", field)
        append = os.path.exists(cal.tables["polleak"])
        cal.casa.polcal(
            vis=cal.vis,
            caltable=cal.tables["polleak"],
            field=field,
            solint="inf",
            combine="scan",
            poltype="Df",
            refant=cal.refant,
            minsnr=2.0,
            minblperant=1,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
        )
    if not os.path.isdir(cal.tables["polleak"]):
        cal.logger.critical("Problem with polarization leakage calibration")
        raise ValueError("Problem with polarization leakage calibration")
    cal.logger.info("Done.")


def polangle_table(cal):
    """
    Derive polarization angle calibration table.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: Nothing
    """
    fields = cal.pol_angle_cals

    cal.logger.info(
        "Calculating the polarization angle calibration table on scan "
        "timescales for polarization angle calibrators ..."
    )
    if os.path.isdir(cal.tables["polangle"]):
        cal.casa.rmtables(cal.tables["polangle"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("polangle", field)
        append = os.path.exists(cal.tables["polangle"])
        cal.casa.polcal(
            vis=cal.vis,
            caltable=cal.tables["polangle"],
            field=field,
            solint="inf",
            combine="scan",
            poltype="Xf",
            refant=cal.refant,
            minsnr=2.0,
            minblperant=1,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
        )
    if not os.path.isdir(cal.tables["polangle"]):
        cal.logger.critical("Problem with polarization angle calibration")
        raise ValueError("Problem with polarization angle calibration")
    cal.logger.info("Done.")
