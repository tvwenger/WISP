"""
caltables.py - Generate calibration tables

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


def set_cal_models(cal):
    """
    Set the flux and polarization calibrator models

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: Nothing
    """
    cal.logger.info("Setting flux and polarization calibrator flux models...")
    for calibrator in list(
        set(cal.flux_cals + cal.pol_leak_cals + cal.pol_angle_cals)
    ):
        # If calibrator model is supplied in config, use that
        cal_models = cal.config.get("Calibrator Models", "Name").splitlines()
        if calibrator in cal_models:
            idx = cal_models.index(calibrator)

            # get reference frequency
            reffreq = cal.config.get(
                "Calibrator Models", "Reference Frequency"
            ).splitlines()[idx]

            # get stokes I flux density
            fluxdensity = cal.config.get(
                "Calibrator Models", "Flux Density"
            ).splitlines()
            fluxdensity = [float(fluxdensity[idx]), 0, 0, 0]

            # get spectral index coefficients
            spix = cal.config.get(
                "Calibrator Models", "Spectral Index Coefficients"
            ).splitlines()[idx]
            spix = [float(i) for i in spix.split(",")]

            # get polarization fraction coefficients
            polindex = cal.config.get(
                "Calibrator Models", "Polarization Fraction Coefficients"
            ).splitlines()[idx]
            polindex = [float(i) for i in polindex.split(",")]

            # get polarization angle coefficients
            polangle = cal.config.get(
                "Calibrator Models", "Polarization Angle Coefficients"
            ).splitlines()[idx]
            polangle = [float(i) for i in polangle.split(",")]

            # Run setjy in manual mode
            cal.casa.setjy(
                vis=cal.vis,
                field=calibrator,
                scalebychan=True,
                standard="manual",
                fluxdensity=fluxdensity,
                spix=spix,
                reffreq=reffreq,
                polindex=polindex,
                polangle=polangle,
            )

        # Otherwise, use CASA model if this is a flux calibrator
        elif calibrator in cal.flux_cals:
            cal.casa.setjy(vis=cal.vis, field=calibrator, scalebychan=True)

        else:
            cal.logger.warn(
                "Not setting a flux model for {0}".format(calibrator)
            )
    cal.logger.info("Done.")


def initial_tables(cal):
    """
    Set the antenna position, gaincurve, and opacity tables.

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: Nothing
    """
    # Correct antenna positions
    if cal.antpos:
        caltable = "antpos.cal"
        cal.logger.info("Calculating antenna position correction...")
        if os.path.isdir(caltable):
            cal.casa.rmtables(caltable)
        cal.casa.gencal(vis=cal.vis, caltable=caltable, caltype="antpos")
        if os.path.isdir(caltable):
            cal.tables["antpos"] = caltable
        cal.logger.info("Done.")

    # Correct for gaincurve and antenna efficiencies
    if cal.gaincurve:
        caltable = "gaincurve.cal"
        cal.logger.info(
            "Calculating gain curve and antenna " "efficiencies..."
        )
        if os.path.isdir(caltable):
            cal.casa.rmtables(caltable)
        cal.casa.gencal(vis=cal.vis, caltable=caltable, caltype="gceff")
        if os.path.isdir(caltable):
            cal.tables["gaincurve"] = caltable
        cal.logger.info("Done.")

    # Correct for atmospheric opacity
    if cal.opacity:
        caltable = "opacity.cal"
        cal.logger.info("Calculating opacity correction...")
        if os.path.isdir(caltable):
            cal.casa.rmtables(caltable)
        myTau = cal.casa.plotweather(vis=cal.vis, doPlot=False)
        allspws = ",".join([str(spw) for spw in range(len(myTau))])
        cal.casa.gencal(
            vis=cal.vis,
            caltable=caltable,
            caltype="opac",
            parameter=myTau,
            spw=allspws,
        )
        if os.path.isdir(caltable):
            cal.tables["opacity"] = caltable
        cal.logger.info("Done.")


def prebandpass_primary_tables(cal, use_smodel=False):
    """
    Derive the delay and phase calibration tables for the primary calibrators
    before bandpass calibration

    Inputs:
        cal :: Calibration object
            The calibration object
        use_smodel :: boolean
            If True, use the Stokes parameter model

    Returns: Nothing
    """
    fields = cal.pri_cals

    # Delay calibration
    cal.logger.info(
        "Calculating delay calibration table for primary calibrators..."
    )
    if os.path.isdir(cal.tables["delays"]):
        cal.casa.rmtables(cal.tables["delays"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("delays", field)
        append = os.path.exists(cal.tables["delays"])
        smodel = None
        if use_smodel:
            smodel = cal.smodels.get(field)
        cal.casa.gaincal(
            vis=cal.vis,
            caltable=cal.tables["delays"],
            field=field,
            refant=cal.refant,
            gaintype="K",
            solint="inf",
            combine="scan",
            minblperant=1,
            parang=cal.calpol,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
            smodel=smodel,
        )
    if not os.path.isdir(cal.tables["delays"]):
        cal.logger.critical("Problem with delay calibration")
        raise ValueError("Problem with delay calibration!")

    # Integration timescale phase calibration
    cal.logger.info(
        "Calculating phase calibration table on integration timescales for "
        "primary calibrators..."
    )
    if os.path.isdir(cal.tables["phase_int0"]):
        cal.casa.rmtables(cal.tables["phase_int0"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("phase_int0", field)
        append = os.path.exists(cal.tables["phase_int0"])
        smodel = None
        if use_smodel:
            smodel = cal.smodels.get(field)
        cal.casa.gaincal(
            vis=cal.vis,
            caltable=cal.tables["phase_int0"],
            field=field,
            solint=cal.solint,
            calmode="p",
            refant=cal.refant,
            gaintype="G",
            minsnr=2.0,
            minblperant=1,
            parang=cal.calpol,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
            smodel=smodel,
        )
    if not os.path.isdir(cal.tables["phase_int0"]):
        cal.logger.critical(
            "Problem with integration-timescale phase calibration"
        )
        raise ValueError(
            "Problem with integration-timescale phase calibration!"
        )
    cal.logger.info("Done.")


def bandpass_table(cal, use_smodel=False):
    """
    Derive bandpass calibration table

    Inputs:
        cal :: Calibration object
            The calibration object
        use_smodel :: boolean
            If True, use the Stokes parameter model

    Returns: Nothing
    """
    fields = cal.pri_cals

    cal.logger.info(
        "Calculating bandpass calibration table for primary calibrators..."
    )
    if os.path.isdir(cal.tables["bandpass"]):
        cal.casa.rmtables(cal.tables["bandpass"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("bandpass", field)
        append = os.path.exists(cal.tables["bandpass"])
        smodel = None
        if use_smodel:
            smodel = cal.smodels.get(field)
        cal.casa.bandpass(
            vis=cal.vis,
            caltable=cal.tables["bandpass"],
            field=field,
            refant=cal.refant,
            fillgaps=1024,
            solint="inf",
            combine="scan",
            solnorm=True,
            minblperant=1,
            parang=cal.calpol,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
            smodel=smodel,
        )
    if not os.path.isdir(cal.tables["bandpass"]):
        cal.logger.critical("Problem with bandpass calibration")
        raise ValueError("Problem with bandpass calibration!")
    cal.logger.info("Done.")


def gain_tables(cal, use_smodel=False):
    """
    Derive phase calibration tables on integration and scan timescales,
    and the amplitude calibration table on scan timescales.

    Inputs:
        cal :: Calibration object
            The calibration object
        use_smodel :: boolean
            If True, use the Stokes parameter model

    Returns: Nothing
    """
    fields = cal.calibrators

    # integration timescale phase calibration
    cal.logger.info(
        "Re-calculating the phase calibration table on integration timescales "
        "for all calibrators..."
    )
    if os.path.isdir(cal.tables["phase_int1"]):
        cal.casa.rmtables(cal.tables["phase_int1"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("phase_int1", field)
        append = os.path.exists(cal.tables["phase_int1"])
        smodel = None
        if use_smodel:
            smodel = cal.smodels.get(field)
        cal.casa.gaincal(
            vis=cal.vis,
            caltable=cal.tables["phase_int1"],
            field=field,
            solint=cal.solint,
            calmode="p",
            refant=cal.refant,
            gaintype="G",
            minsnr=2.0,
            minblperant=1,
            parang=cal.calpol,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
            smodel=smodel,
        )
    if not os.path.isdir(cal.tables["phase_int1"]):
        cal.logger.critical(
            "Problem with integration-timescale phase calibration"
        )
        raise ValueError(
            "Problem with integration-timescale phase calibration!"
        )
    cal.logger.info("Done.")

    # scan timescale phase calibration
    cal.logger.info(
        "Calculating the phase calibration table on "
        "scan timescales for all calibrators..."
    )
    if os.path.isdir(cal.tables["phase_scan"]):
        cal.casa.rmtables(cal.tables["phase_scan"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("phase_scan", field)
        append = os.path.exists(cal.tables["phase_scan"])
        smodel = None
        if use_smodel:
            smodel = cal.smodels.get(field)
        cal.casa.gaincal(
            vis=cal.vis,
            caltable=cal.tables["phase_scan"],
            field=field,
            solint="inf",
            calmode="p",
            refant=cal.refant,
            gaintype="G",
            minsnr=2.0,
            minblperant=1,
            parang=cal.calpol,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
            smodel=smodel,
        )
    if not os.path.isdir(cal.tables["phase_scan"]):
        cal.logger.critical("Problem with scan-timescale phase calibration")
        raise ValueError("Problem with scan-timescale phase calibration!")
    cal.logger.info("Done.")

    # scan timescale amplitude corrections
    cal.logger.info(
        "Calculating the amplitude calibration table "
        "on scan timescales for all calibrators..."
    )
    if os.path.isdir(cal.tables["amplitude"]):
        cal.casa.rmtables(cal.tables["amplitude"])
    for field in fields:
        gaintables, gainfields, spwmaps = cal.gaintables("amplitude", field)
        append = os.path.exists(cal.tables["amplitude"])
        smodel = None
        if use_smodel:
            smodel = cal.smodels.get(field)
        cal.casa.gaincal(
            vis=cal.vis,
            caltable=cal.tables["amplitude"],
            field=field,
            solint="inf",
            calmode="ap",
            refant=cal.refant,
            minsnr=2.0,
            minblperant=1,
            parang=cal.calpol,
            gaintable=gaintables,
            gainfield=gainfields,
            spwmap=spwmaps,
            append=append,
            smodel=smodel,
        )
    if not os.path.isdir(cal.tables["amplitude"]):
        cal.logger.critical("Problem with amplitude calibration")
        raise ValueError("Problem with amplitude calibration!")
    cal.logger.info("Done.")


def flux_table(cal):
    """
    Derive incremental flux scale table

    Inputs:
        cal :: Calibration object
            The calibration object

    Returns: Nothing
    """
    field = ",".join(cal.flux_cals)
    cal.logger.info("Calculating the flux scale calibration table...")
    if os.path.isdir(cal.tables["flux"]):
        cal.casa.rmtables(cal.tables["flux"])
    cal.casa.fluxscale(
        vis=cal.vis,
        caltable=cal.tables["amplitude"],
        fluxtable=cal.tables["flux"],
        reference=field,
        incremental=True,
    )
    if not os.path.isdir(cal.tables["flux"]):
        cal.logger.critical("Problem with flux scale alibration")
        raise ValueError("Problem with flux scale calibration!")
    cal.logger.info("Done.")
