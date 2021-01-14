"""
flagging.py - WISP Flagging Routines

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

import time
import pickle

import numpy as np


def preliminary_flagging(cal, casa):
    """
    Perform preliminary flagging: shadowed antennas, quack,
    flags from configuration file, interpolations from
    configuration file, then extend flags as necessary.

    Inputs:
        cal :: Calibration object
            The calibration object
        casa :: CASA namespace
            CASA namespace

    Returns: Nothing
    """
    cal.save_flags("starting_flags")

    # Flag shadowed antennas
    cal.logger.info("Flagging shadowed antennas...")
    casa.flagdata(
        vis=cal.vis,
        mode="shadow",
        tolerance=cal.shadow_tolerance,
        flagbackup=False,
        extendflags=False,
    )
    cal.logger.info("Done.")

    # Flag the beginning of each scan
    if cal.quack_interval > 0:
        cal.logger.info(
            "Flagging the first %f of each scan...", cal.quack_interval
        )
        casa.flagdata(
            vis=cal.vis,
            mode="quack",
            quackinterval=cal.quack_interval,
            flagbackup=False,
            extendflags=False,
        )
        cal.logger.info("Done.")

    # Flag scans from configuration file
    scan = cal.config.get("Flags", "Scan")
    if scan != "":
        cal.logger.info("Flagging scans from configuration file: " "%s", scan)
        casa.flagdata(
            vis=cal.vis,
            mode="manual",
            scan=scan,
            flagbackup=False,
            extendflags=False,
        )
        cal.logger.info("Done.")

    # Flag antennas from configuration file
    antenna = cal.config.get("Flags", "Antenna")
    if antenna != "":
        cal.logger.info(
            "Flagging antennas from configuration " "file: %s", antenna
        )
        casa.flagdata(
            vis=cal.vis,
            mode="manual",
            antenna=antenna,
            flagbackup=False,
            extendflags=False,
        )
        cal.logger.info("Done.")

    # Flag spectral windows from configuration file
    spws = cal.config.get("Flags", "Spectral Window")
    if spws != "":
        cal.logger.info(
            "Flagging spectral windows from " "configuration file: %s",
            spws,
        )
        casa.flagdata(
            vis=cal.vis,
            mode="manual",
            spw=spws,
            flagbackup=False,
            extendflags=False,
        )
        cal.logger.info("Done.")

    # Flag line channels from configuration file
    badchans = cal.config.get("Flags", "Line Channels")
    badchans = badchans.split(",")
    if badchans[0] != "":
        cal.logger.info(
            "Flagging line channels from configuration " "file: %s",
            ";".join(badchans),
        )
        badchans = ";".join(badchans)
        line_spws = ",".join(
            [i + ":" + badchans for i in cal.line_spws.split(",")]
        )
        casa.flagdata(
            vis=cal.vis,
            mode="manual",
            spw=line_spws,
            flagbackup=False,
            extendflags=False,
        )
        cal.logger.info("Done.")

    # Flag continuum channels from configuration file
    badchans = cal.config.get("Flags", "Continuum Channels")
    badchans = badchans.split(",")
    if badchans[0] != "":
        cal.logger.info(
            "Flagging continuum channels from " "configuration file: %s",
            ";".join(badchans),
        )
        badchans = ";".join(badchans)
        cont_spws = ",".join(
            [i + ":" + badchans for i in cal.cont_spws.split(",")]
        )
        casa.flagdata(
            vis=cal.vis,
            mode="manual",
            spw=cont_spws,
            flagbackup=False,
            extendflags=False,
        )
        cal.logger.info("Done.")

    # Interpolate through bad line channels
    badchans = cal.config.get("Interpolate", "Line Channels")
    badchans = badchans.split(",")
    if badchans[0] != "":
        cal.logger.info(
            "Interpolating through line channels from "
            "configuration file: %s",
            ";".join(badchans),
        )
        line_spws = [int(i) for i in cal.line_spws.split(",")]
        badchans = np.array([int(i) for i in badchans])
        cal.interpolate_channels(line_spws, badchans)
        cal.logger.info("Done.")

    # Interpolate through bad continuum channels
    badchans = cal.config.get("Interpolate", "Continuum Channels")
    badchans = badchans.split(",")
    if badchans[0] != "":
        cal.logger.info(
            "Interpolating through continuum channels "
            "from configuration file: %s",
            ";".join(badchans),
        )
        cont_spws = [int(i) for i in cal.cont_spws.split(",")]
        badchans = np.array([int(i) for i in badchans])
        cal.interpolate_channels(cont_spws, badchans)
        cal.logger.info("Done")

    # Extend the flags
    cal.logger.info("Extending flags...")
    casa.flagdata(
        vis=cal.vis,
        mode="extend",
        extendpols=True,
        growtime=90.0,
        growfreq=90.0,
        growaround=True,
        flagbackup=False,
    )
    cal.logger.info("Done.")
    cal.save_flags("preliminary")


def auto_flag(cal, casa, fields, extend=False):
    """
    Perform automatic flagging using the rflag algorithm on
    calibrated visiblities, then optionally extend the flags.

    Inputs:
        cal :: Calibration object
            The calibration object
        casa :: CASA namespace
            CASA namespace
        field :: list of strings
            The field(s) to flag
        extend :: boolean
            If True, extend the flags following automatic flagging

        Returns: Nothing
    """
    cal.logger.info("Running rflag on corrected data column...")
    casa.flagdata(
        vis=cal.vis,
        mode="rflag",
        field=",".join(fields),
        flagbackup=False,
        datacolumn="corrected",
        extendflags=False,
    )
    cal.logger.info("Done.")

    # Extend the flags
    if extend:
        cal.logger.info("Extending flags...")
        casa.flagdata(
            vis=cal.vis,
            mode="extend",
            extendpols=True,
            growtime=90.0,
            growfreq=90.0,
            growaround=True,
            flagbackup=False,
        )
        cal.logger.info("Done.")

    # Save the flags
    cal.save_flags("autoflag")


def manual_flag(cal, fieldtype):
    """
    Interactively plot and flag.

    Inputs:
        cal :: Calibration object
            The calibration object
        fieldtype :: string
            Either 'calibrator' or 'science'

    Returns: Nothing
    """
    fname = "science_plots.pkl"
    datacolumn = "corrected"
    if fieldtype == "calibrator":
        fname = "calibrator_plots.pkl"
        # get datacolumn
        datacolumn = input("Datacolumn? (data or corrected) ")

    # Read plot list from pickle object
    cal.logger.info("Reading plot list from pickle...")
    with open(fname, "r") as fin:
        plots = pickle.load(fin)
    num_plots = len(plots)
    cal.logger.info("Done.")
    #
    # Prompt user with menu
    #
    cal.logger.info("Please inspect PDF and then perform manual flagging.")
    while True:
        print("f :: flag some data")
        print("<number> :: open interactive plot with this id")
        print("quit :: end this flagging session")
        answer = input()

        # Flag some data
        if answer.lower() == "f":
            flag(cal, fieldtype)

        # Stop flagging
        elif answer.lower() == "quit":
            break

        # Generate plotms figure
        else:
            try:
                plotid = int(answer)
            except ValueError:
                print("Invalid Plot ID")
                continue
            if plotid >= num_plots:
                print("Invalid Plot ID")
                continue
            title = "PlotID: {0} Field: {1}".format(
                plotid, plots[plotid]["field"]
            )
            cal.casa.plotms(
                vis=cal.vis,
                xaxis=plots[plotid]["xaxis"],
                yaxis=plots[plotid]["yaxis"],
                field=plots[plotid]["field"],
                ydatacolumn=datacolumn,
                iteraxis="spw",
                title=title,
                coloraxis="baseline",
                correlation=cal.correlation,
                avgchannel=plots[plotid]["avgchannel"],
                avgtime=plots[plotid]["avgtime"],
            )

    # Save the flags
    cal.save_flags("manualflag")


def flag(cal, fieldtype):
    """
    Interactively flag some data.

    Inputs:
        cal :: Calibration object
            The calibration object
        fieldtype :: string
            Either 'calibrator' or 'science'

    Returns: Nothing
    """
    fields = cal.sci_targets
    if fieldtype == "calibrator":
        fields = cal.all_fields

    # Build list of flag commands
    flag_commands = []
    while True:

        # Prompt user for attributes to flag
        print("Field? Empty = {0}".format(",".join(fields)))
        field = input()
        if field == "":
            field = ",".join(fields)
        if field not in ",".join(fields):
            print("{0} is not in {1}".format(field, fields))
            continue
        print(
            "Scan? Empty = all scans (ex. 0 to flag scan 0, 1~3 "
            "to flag scans 1, 2, and 3)"
        )
        scan = input()
        print(
            "Spectral window and channels? Empty = all spws "
            "(ex. 2:100~150;200:250,5:1200~1210 is spw 2, chans "
            "100 to 150 and 200 to 250, and spw 5 chans 1200 to "
            "1210)"
        )
        spw = input()
        print("Time range? Empty = all times (ex. " "10:23:45~10:23:55)")
        timerange = input()
        print(
            "Antenna or baseline? Empty = all antennas/baselines "
            "(ex. CA01 to flag ant 1 or CA01&CA02 to flag 1-2 "
            "baseline)"
        )
        antenna = input()
        print("Correlation? Empty = all correlations (i.e. XX,YY)")
        correlation = input()

        # Build flag command
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
        flag_commands.append(" ".join(parts))

        # Confirm with user, or append more flag commands
        print("Will execute:")
        print(
            "flagdata(vis='{0}',mode='list',flagbackup=False,"
            "extendflags=False,".format(cal.vis)
        )
        for icmd, cmd in enumerate(flag_commands):
            if len(flag_commands) == 1:
                print('         inpfile=["{0}"])'.format(cmd))
            elif icmd == 0:
                print('         inpfile=["{0}",'.format(cmd))
            elif icmd == len(flag_commands) - 1:
                print('                  "{0}"])'.format(cmd))
            else:
                print('                  "{0}",'.format(cmd))
        print("Proceed [y/n] or add another flag command [a]?")
        answer = input()

        # Execute flag command
        if answer.lower() == "y":
            cal.logger.info("Executing:")
            cal.logger.info(
                "flagdata(vis='%s',mode='list',"
                "flagbackup=False,extendflags=False,",
                cal.vis,
            )
            for icmd, cmd in enumerate(flag_commands):
                if len(flag_commands) == 1:
                    cal.logger.info('         inpfile=["%s"])', cmd)
                elif icmd == 0:
                    cal.logger.info('         inpfile=["%s",', cmd)
                elif icmd == len(flag_commands) - 1:
                    cal.logger.info('                  "%s"])', cmd)
                else:
                    cal.logger.info('                  "%s",', cmd)

            # Save flag command to manual flags list
            with open("manual_flags.txt", "a") as fout:
                cur_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
                for cmd in flag_commands:
                    fout.write("{0}: {1}".format(cur_time, cmd) + "\n")

            # Execute
            cal.casa.flagdata(
                vis=cal.vis,
                mode="list",
                flagbackup=False,
                extendflags=False,
                inpfile=flag_commands,
            )
            break

        # Append another flag command
        elif answer.lower() == "a":
            continue

        # Quit flagging
        else:
            print("Aborting...")
            break
