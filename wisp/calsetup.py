"""
calsetup.py - Initialize a WISP calibration object

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

import numpy as np
from .utils import natural_sort


def get_calibrators(cal, caltype):
    """
    Get calibrators from the configuration file or the listobs file.

    Inputs:
        cal :: Calibration object
            The calibration object
        caltype :: string
            The calibrator type. One of: "Primary", "Secondary", "Flux",
            "PolLeakage", or "PolAngle"

    Returns: cals
        cals :: list of strings
            Calibrator field names
    """
    # First try getting calibrators
    cals = cal.config.get(
        "Calibrators", "{0} Calibrators".format(caltype)
    ).splitlines()
    if cals:
        return [field for field in cals if field in cal.all_fields]

    # Get calibrators from listobs
    cals = []

    # map calibrator type to expected listobs item
    calmap = {
        "Primary": "CALIBRATE_BANDPASS",
        "Secondary": "CALIBRATE_PHASE",
        "Flux": "CALIBRATE_FLUX",
        "PolLeakage": "CALIBRATE_POL_LEAKAGE",
        "PolAngle": "CALIBRATE_POL_ANGLE",
    }
    with open("listobs.txt", "r") as fin:
        for line in fin:
            for field in cal.all_fields:
                if field in line and calmap[caltype] in line:
                    cals.append(field)
    return list(set(cals))


def get_scan_fields(cal):
    """
    Get field names for each scan.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: fields
        fields :: list of strings
            Field names associated with each scan, in scan order.
    """
    # get scan numbers
    cal.casa.ms.open(cal.vis)
    scans = cal.casa.ms.getscansummary()
    cal.casa.ms.close()
    scan_nums = natural_sort(scans.keys())

    # get field names
    cal.casa.msmd.open(cal.vis)
    fieldnames = cal.casa.msmd.fieldnames()
    cal.casa.msmd.close()
    return [fieldnames[scans[num]["0"]["FieldId"]] for num in scan_nums]


def assign_secondary_calibrators(cal):
    """
    Determine which secondary calibrator to use for each science target,
    based on which calibrator is closest in time.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: assignments
        assignments :: dictionary
            For each science target field name, the assigned secondary
            calibrator field name.
    """
    # Get field names associated with each scan
    scan_fields = get_scan_fields(cal)

    # get indicies of secondary calibrators
    cal_inds = np.array(
        [i for i in range(len(scan_fields)) if scan_fields[i] in cal.sec_cals]
    )

    # Assign closest secondary calibrator (in time) to each science target
    science_calibrators = {}
    cal.logger.info("{0:15} {1:15}".format("Target", "Calibrator"))
    for sci_field in scan_fields:
        # skip calibrators and those already assigned
        if sci_field in cal.calibrators or sci_field in science_calibrators:
            continue

        # Get indicies of field in all scans
        sci_inds = np.array(
            [i for i in range(len(scan_fields)) if scan_fields[i] == sci_field]
        )

        calib = None
        for ind in sci_inds:
            # get closest before calibrator
            try:
                closest_before = scan_fields[
                    cal_inds[cal_inds < ind][
                        np.argmin(np.abs(cal_inds[cal_inds < ind] - ind))
                    ]
                ]
            except ValueError:
                # none before
                closest_before = None

            # get closest after calibrator
            try:
                closest_after = scan_fields[
                    cal_inds[cal_inds > ind][
                        np.argmin(np.abs(cal_inds[cal_inds > ind] - ind))
                    ]
                ]
            except ValueError:
                # none after
                closest_after = None

            # If we've already found a calibrator for this source,
            # check that we found the same calibrator again.
            if calib is not None:
                if (closest_before == closest_after) and (
                    closest_before == calib
                ):
                    # good
                    continue
                elif (closest_before == closest_after) and (
                    closest_before != calib
                ):
                    # update
                    calib = closest_before
                elif (
                    (closest_before is not None) and (closest_before == calib)
                ) or (
                    (closest_after is not None) and (closest_after == calib)
                ):
                    # good
                    continue
                else:
                    raise ValueError(
                        "%s has different closest cals: "
                        "calib: %s before: %s after: %s",
                        sci_field,
                        calib,
                        closest_before,
                        closest_after,
                    )
            else:
                if closest_before is None and closest_after is None:
                    raise ValueError(
                        "No good calibrator {0}".format(sci_field)
                    )
                elif closest_before == closest_after:
                    # good
                    calib = closest_before
                elif closest_before is not None:
                    calib = closest_before
                else:
                    calib = closest_after
        science_calibrators[sci_field] = calib
        cal.logger.info("{0:15} {1:15}".format(sci_field, calib))
    return science_calibrators
