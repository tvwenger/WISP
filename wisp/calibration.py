"""
calibration.py - WISP Calibration Pipeline Object

Calibrate a measurement set by scripting CASA tasks and generating
diagnostic plots.

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
Trey V. Wenger November 2018 - V1.0

Trey V. Wenger December 2018 - V1.1
    Removed TFCROP - it was causing problems with VLA data.
    Added preliminary scan flagging, extend flags, and tunable
    parameters for shadow tolerance and quack interval.
    Added opacity, gaincurve, and antenna position calibrations.
    Changed plotcal figures to use plotms for generation.

Trey V. Wenger August 2019 - V2.0
    Added polarization leakage calibration.
    Re-designed to OOP framework.
    Speed up interpolation.

Trey V. Wenger - April 2021 - v3.0
    Improve code readability.
    Restructure polarization calibration.
    Add support for polarization angle calibration.

Trey V. Wenger - April 2021 - v3.1
    Fix bug in calibration table order and applycal parang
"""

import os
import time

from .calsetup import (
    get_calibrators,
    get_scan_fields,
    get_spw_corrs,
    assign_secondary_calibrators,
)
from .caltables import (
    flux_table,
    set_cal_models,
    initial_tables,
    prebandpass_primary_tables,
    bandpass_table,
    gain_tables,
)

from .poltables import (
    crosshand_delays_table,
    polangle_table,
    polleak_table,
)

from .utils import get_smodels, natural_sort


class Calibration:
    """
    The Calibration object handles the calibration steps for a
    measurement set.
    """

    def __init__(
        self,
        casa,
        vis,
        logger,
        config,
        refant,
        shadow_tolerance=0.0,
        quack_interval=10.0,
        antpos=True,
        gaincurve=True,
        opacity=True,
        calpol=False,
        leakage_solint="inf",
        calwt=True,
        solint="int",
    ):
        """
        Create a new Calibration object. Get reference antenna,
        create listobs file, find line and continuum spectral windows,
        generate lists of calibrators, and create needed directories.

        Inputs:
            casa :: CASA namespace
                CASA namespace
            vis :: string
                The masurement set
            logger :: logging.Logger object
                The logging object we're using
            config :: config.ConfigParser object
                The config parser we're using
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
                if True, compute antenna position corrections
                (only for VLA)
            gaincurve :: boolean
                if True, compute gain curve and antenna efficiency
                corrections (only for VLA)
            opacity :: boolean
                if True, compute opacity corrections
            calpol :: boolean
                if True, calibrate polarization
            leakge_solint :: string
                Solution interval for D-term calibration (e.g., solint="inf,2MHz")
            calwt :: boolean
                if True, apply calibration weights to data
            solint :: string
                The solution interval for short timescale phase corrections.
                Default is 'int' for integration. (e.g., solint='10s')

        Returns: calibration
            calibration :: calibration.Calibration object
                a new Calibration object
        """
        self.casa = casa
        self.vis = vis
        self.logger = logger
        self.config = config
        self.refant = refant
        self.shadow_tolerance = shadow_tolerance
        self.quack_interval = quack_interval
        self.antpos = antpos
        self.gaincurve = gaincurve
        self.opacity = opacity
        self.calpol = calpol
        self.leakage_solint = leakage_solint
        self.calwt = calwt
        self.solint = solint
        self.logger.info("Initializing Calibration object.")

        # storage for flux models returned by fluxscale
        self.flux_models = None

        # Initialize calibration tables
        self.tables = {
            "antpos": "antpos.cal",
            "gaincurve": "gaincurve.cal",
            "opacity": "opacity.cal",
            "delays": "delays.Kcal",
            "phase_int0": "phase_int0.Gcal",
            "bandpass": "bandpass.Bcal",
            "phase_int1": "phase_int1.Gcal",
            "phase_scan": "phase_scan.Gcal",
            "amplitude": "amplitude.Gcal",
            "flux": "flux.cal",
            "crosshand_delays": "crosshand_delays.Kcal",
            "polleak": "polleak.Dcal",
            "polangle": "polangle.Xcal",
        }

        # Generate listobs file
        listfile = "listobs.txt"
        if not os.path.isfile(listfile):
            self.logger.info("Generating listobs file...")
            self.casa.listobs(vis=self.vis, listfile=listfile)
            self.logger.info("Done.")

        # Get continuum and line spws from configuration file
        self.cont_spws = (
            self.config.get("Spectral Windows", "Continuum").strip().split(",")
        )
        self.cont_spws = [spw.strip() for spw in self.cont_spws if spw.strip() != ""]
        self.line_spws = self.config.get("Spectral Windows", "Line").strip().split(",")
        self.line_spws = [spw.strip() for spw in self.line_spws if spw.strip() != ""]
        self.logger.info("Found continuum spws: %s", ",".join(self.cont_spws))
        self.logger.info("Found line spws: %s", ",".join(self.line_spws))

        # Get feed orientation for each spectral window
        self.spw_corrs = get_spw_corrs(self)
        self.orientation = None
        self.logger.info("Spectral window observed correlations:")
        for spw in natural_sort(self.spw_corrs.keys()):
            self.logger.info("%s: %s", spw, self.spw_corrs[spw])
            if "RR" in self.spw_corrs[spw] or "LL" in self.spw_corrs[spw]:
                orientation = "circular"
            elif "XX" in self.spw_corrs[spw] or "YY" in self.spw_corrs[spw]:
                orientation = "linear"
            else:
                raise ValueError("Unknown correlation: {0}".format(self.spw_corrs[spw]))
            if self.orientation is None:
                self.orientation = orientation
            elif self.orientation != orientation:
                raise ValueError("Correlation mismatch!")

        # get spectral windows with and without cross correlations
        self.all_spws = []
        self.corr_spws = []
        self.nocorr_spws = []
        for i in natural_sort(self.spw_corrs.keys()):
            self.all_spws.append(i)
            if "RL" in self.spw_corrs[i] or "XY" in self.spw_corrs[i]:
                self.corr_spws.append(i)
            else:
                self.nocorr_spws.append(i)
        if self.calpol:
            self.logger.info(
                "Polarization calibration will be applied to these "
                "spectral windows: %s",
                ",".join(self.corr_spws),
            )

        # get unique field names
        self.logger.info("Looking for field names...")
        self.all_fields = list(set(get_scan_fields(self)))
        self.logger.info("Found fields: %s", ", ".join(self.all_fields))

        # Get primary calibrator fields
        self.pri_cals = get_calibrators(self, "Primary")
        self.logger.info("Primary calibrators: %s", ", ".join(self.pri_cals))

        # Get Secondary calibrator fields
        self.sec_cals = get_calibrators(self, "Secondary")
        self.logger.info("Secondary calibrators: %s", ", ".join(self.sec_cals))

        # Get flux calibrator fields
        self.flux_cals = get_calibrators(self, "Flux")
        self.logger.info("Flux calibrators: %s", ", ".join(self.flux_cals))

        # Get polarization leakage calibrators
        self.pol_leak_cals = get_calibrators(self, "PolLeakage")

        # Get polarization angle calibrators
        self.pol_angle_cals = get_calibrators(self, "PolAngle")

        # Check that polarization calibrators exist if needed
        if self.calpol:
            if not self.pol_leak_cals:
                raise ValueError("No polarization leakge calibrators found")
            self.logger.info(
                "Polarization leakage calibrators: %s",
                ", ".join(self.pol_leak_cals),
            )
            if self.orientation == "circular":
                if not self.pol_angle_cals:
                    raise ValueError("No polarization angle calibrators found.")
                self.logger.info(
                    "Polarization angle calibrators: %s",
                    ", ".join(self.pol_angle_cals),
                )

        # Compile calibrator fields
        self.calibrators = list(
            set(
                self.pri_cals
                + self.sec_cals
                + self.flux_cals
                + self.pol_leak_cals
                + self.pol_angle_cals
            )
        )

        # Get science targets
        self.sci_targets = [
            field for field in self.all_fields if field not in self.calibrators
        ]
        self.logger.info("Science targets: %s", ", ".join(self.sci_targets))

        # Determine which secondary calibrator to user for each
        # science target
        self.logger.info("Identifying secondary calibrators for each science target...")
        self.science_calibrators = assign_secondary_calibrators(self)

        # create directories for figures
        if not os.path.isdir("calibrator_figures"):
            os.makedirs("calibrator_figures")
        if not os.path.isdir("science_figures"):
            os.makedirs("science_figures")
        if not os.path.isdir("plotcal_figures"):
            os.makedirs("plotcal_figures")

    def save_flags(self, label):
        """
        Save the flag state with a given label and the current
        time.

        Inputs:
          label :: string
            A label to add to the version name

        Returns: Nothing
        """
        self.logger.info("Saving flag state...")
        cur_time = time.strftime("%Y %m %d %H %M %S", time.gmtime())
        versionname = "{0} {1}".format(label, cur_time)
        self.casa.flagmanager(vis=self.vis, mode="save", versionname=versionname)
        self.logger.info("Done")

    def gaintables(self, step, field, poltables=True):
        """
        Get the calibration tables up to a given calibration step, excluding
        missing optional tables. Also return the associated fields and spectral
        window mappings.
            Available steps:
                delays, phase_int0, bandpass, phase_int1, phase_scan,
                amplitude, crosshand_delays, polleak, polangle, apply

        Inputs:
            step :: string
                The step for which to get the preceeding calibration tables.
                One of the steps from the list above.
            field :: string
                The field name to which these calibration tables will be
                applied. This parameter sets gainfields for the phase and
                amplitude tables to be either <field> if <field> is a
                calibrator or the associated calibrator <field> if <field>
                is a science target.
            poltables :: boolean
                If True, include polarization tables if self.calpol is True

        Returns: gaintables, gainfields, spwmaps
            gaintables :: list of strings
                The calibration tables required for this step.
            gainfields :: list of strings
                The field names to select from each calibration table.
            spwmaps :: list of list of integers
                The mapping of spectral windows for each calibration table.
        """
        gaintables = []
        gainfields = []
        spwmaps = []

        # optional tables
        if self.antpos and os.path.exists(self.tables["antpos"]):
            gaintables.append(self.tables["antpos"])
            gainfields.append("")
            spwmaps.append([])
        if self.gaincurve and os.path.exists(self.tables["gaincurve"]):
            gaintables.append(self.tables["gaincurve"])
            gainfields.append("")
            spwmaps.append([])
        if self.opacity and os.path.exists(self.tables["opacity"]):
            gaintables.append(self.tables["opacity"])
            gainfields.append("")
            spwmaps.append([])
        if step == "delays":
            return gaintables, gainfields, spwmaps

        # add delays
        gaintables.append(self.tables["delays"])
        gainfields.append("")
        spwmaps.append([])
        if step == "phase_int0":
            return gaintables, gainfields, spwmaps

        # add phase_int0 for bandpass
        if step == "bandpass":
            return (
                gaintables + [self.tables["phase_int0"]],
                gainfields + [field],
                spwmaps + [[]],
            )

        # add bandpass
        gaintables.append(self.tables["bandpass"])
        gainfields.append("")
        spwmaps.append([])
        if step == "phase_int1" or step == "phase_scan":
            return gaintables, gainfields, spwmaps

        # add phase_int1 or phase_scan
        if field in self.sci_targets:
            gaintables.append(self.tables["phase_scan"])
            gainfields.append(self.science_calibrators[field])
            spwmaps.append([])
        else:
            gaintables.append(self.tables["phase_int1"])
            gainfields.append(field)
            spwmaps.append([])
        if step == "amplitude":
            return gaintables, gainfields, spwmaps

        # add amplitude
        gaintables.append(self.tables["amplitude"])
        if field in self.sci_targets:
            gainfields.append(self.science_calibrators[field])
        else:
            gainfields.append(field)
        spwmaps.append([])
        if step == "crosshand_delays":
            return gaintables, gainfields, spwmaps

        # add polarization calibration tables if necessary
        if self.calpol and poltables:
            # add cross-hand delays
            if os.path.exists(self.tables["crosshand_delays"]):
                gaintables.append(self.tables["crosshand_delays"])
                gainfields.append("")
                spwmaps.append([])
            if step == "polleak":
                return gaintables, gainfields, spwmaps

            # add polarization leakage
            if os.path.exists(self.tables["polleak"]):
                gaintables.append(self.tables["polleak"])
                gainfields.append("")
                spwmaps.append([])
            if step == "polangle":
                return gaintables, gainfields, spwmaps

            # add polarization angle
            if os.path.exists(self.tables["polangle"]):
                gaintables.append(self.tables["polangle"])
                gainfields.append("")
                spwmaps.append([])

        # add flux
        gaintables.append(self.tables["flux"])
        if field in self.sci_targets:
            gainfields.append(self.science_calibrators[field])
        else:
            gainfields.append(field)
        spwmaps.append([])
        if step == "apply":
            return gaintables, gainfields, spwmaps

        # we should never get here
        raise ValueError("Invalid step: {0}".format(step))

    def split_fields(self):
        """
        Split calibrated fields into measurement sets with naming:
        {field_name}_calibrated.ms

        Inputs: Nothing

        Returns: Nothing
        """
        for field in self.all_fields:
            outputvis = "{0}_calibrated.ms".format(field)
            self.logger.info("Splitting %s to %s", field, outputvis)
            self.casa.split(
                vis=self.vis, outputvis=outputvis, field=field, keepflags=False
            )
        self.logger.info("Done!")


def generate_tables(cal):
    """
    Generate the calibration tables.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: Nothing
    """
    # set the model for the flux calibrators
    set_cal_models(cal)

    # generate the antenna position, gaincurve, and opacity tables
    initial_tables(cal)

    # pre-bandpass delay calibration
    prebandpass_primary_tables(cal, use_smodel=False)

    # bandpass
    bandpass_table(cal, use_smodel=False)

    # complex gain
    gain_tables(cal, use_smodel=False)

    # polarization calibration
    if cal.calpol:
        # cross-hand delay calibration
        if cal.orientation == "circular":
            crosshand_delays_table(cal)

        # polarization leakage
        polleak_table(cal)

        # polarization angle
        if cal.orientation == "circular":
            polangle_table(cal)

        if cal.orientation == "linear":
            # apply calibration solutions
            apply_calibration(cal, "calibrator")

            # Estimate polarization from corrected data
            cal.smodels = get_smodels(cal)

            # Re-do calibration using estimated polarization
            prebandpass_primary_tables(cal, use_smodel=True)
            bandpass_table(cal, use_smodel=True)
            gain_tables(cal, use_smodel=True)

    # set the flux scale
    flux_table(cal)


def apply_calibration(cal, fieldtype):
    """
    Apply calibration solutions to either calibrator fields or science targets.

    Inputs:
        cal :: Calibration object
            The calibration object
        fieldtype :: string
            Either 'calibrator' or 'science'

    Returns: Nothing
    """
    if fieldtype == "calibrator":
        fields = cal.calibrators
    elif fieldtype == "science":
        fields = cal.sci_targets
    else:
        raise ValueError("Invalid fieldtype: {0}".format(fieldtype))
    for field in fields:
        if cal.calpol:
            # no polarization spectral windows
            if len(cal.nocorr_spws) > 0:
                gaintables, gainfields, spwmaps = cal.gaintables(
                    "apply", field, poltables=False
                )
                cal.casa.applycal(
                    vis=cal.vis,
                    field=field,
                    spw=",".join(cal.nocorr_spws),
                    calwt=cal.calwt,
                    gaintable=gaintables,
                    gainfield=gainfields,
                    spwmap=spwmaps,
                    parang=cal.calpol,
                    flagbackup=False,
                )

            # polarization spectral windows
            if len(cal.corr_spws) > 0:
                gaintables, gainfields, spwmaps = cal.gaintables(
                    "apply", field, poltables=True
                )
                cal.casa.applycal(
                    vis=cal.vis,
                    field=field,
                    spw=",".join(cal.corr_spws),
                    calwt=cal.calwt,
                    gaintable=gaintables,
                    gainfield=gainfields,
                    spwmap=spwmaps,
                    parang=cal.calpol,
                    flagbackup=False,
                )
        else:
            # all spectral windows
            gaintables, gainfields, spwmaps = cal.gaintables("apply", field)
            cal.casa.applycal(
                vis=cal.vis,
                field=field,
                calwt=cal.calwt,
                gaintable=gaintables,
                gainfield=gainfields,
                spwmap=spwmaps,
                parang=cal.calpol,
                flagbackup=False,
            )
    cal.save_flags("calibrate")
