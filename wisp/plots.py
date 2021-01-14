"""
plots.py - WISP plotting utilities

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
import glob
import pickle

from .utils import generate_pdf


def plotcal_plots(cal, casa):
    """
    Generate diagnostic plots for calibration tables.

    Inputs:
        cal :: Calibration object
            The calibration object
        casa :: CASA namespace
            CASA namespace

    Returns: Nothing
    """

    # Remove previous figures
    fnames = glob.glob("plotcal_figures/*.png")
    for fname in fnames:
        os.remove(fname)

    # Bandpass
    field = ",".join(cal.pri_cals)
    casa.plotms(
        vis=cal.bandpass,
        xaxis="channel",
        yaxis="amplitude",
        field=field,
        iteraxis="spw",
        coloraxis="antenna1",
        title=cal.bandpass.replace("_", r"\_"),
        plotfile="plotcal_figures/0_bandpass.png",
        overwrite=True,
        showgui=False,
        exprange="all",
    )

    # Integration phase
    field = list(
        set(
            cal.pri_cals
            + cal.sec_cals
            + cal.pol_leak_cals
            + cal.pol_angle_cals
        )
    )
    casa.plotms(
        vis=cal.phase_int,
        xaxis="time",
        yaxis="phase",
        field=field,
        iteraxis="spw",
        coloraxis="antenna1",
        title=cal.phase_int.replace("_", r"\_"),
        plotfile="plotcal_figures/1_phase_int.png",
        overwrite=True,
        showgui=False,
        exprange="all",
    )

    # Scan phase
    casa.plotms(
        vis=cal.phase_scan,
        xaxis="time",
        yaxis="phase",
        field=field,
        iteraxis="spw",
        coloraxis="antenna1",
        title=cal.phase_scan.replace("_", r"\_"),
        plotfile="plotcal_figures/2_phase_scan.png",
        overwrite=True,
        showgui=False,
        exprange="all",
    )

    # Amplitude scan
    casa.plotms(
        vis=cal.apcal_scan,
        xaxis="time",
        yaxis="amplitude",
        field=field,
        iteraxis="spw",
        coloraxis="antenna1",
        title=cal.apcal_scan.replace("_", r"\_"),
        plotfile="plotcal_figures/3_apcal_scan.png",
        overwrite=True,
        showgui=False,
        exprange="all",
    )

    # Generate PDF of plotcal figures
    cal.logger.info("Generating PDF...")
    generate_pdf("plotcal_figures")
    cal.logger.info("Done.")


def visibility_plots(cal, casa, fieldtype):
    """
    Generate diagnostic visiblity plots for calibrators.

    Inputs:
        cal :: Calibration object
            The calibration object
        casa :: CASA namespace
            CASA namespace
        fieldtype :: string
            Either 'calibrator' or 'science'

    Returns: Nothing
    """
    # get datacolumn and pickle file
    datacolumn = "corrected"
    location = "science_figures"
    fields = cal.sci_targets
    if fieldtype == "calibrator":
        datacolumn = input("Datacolumn? (data or corrected) ")
        location = "calibrator_figures"
        fields = list(
            set(
                cal.pri_cals
                + cal.sec_cals
                + cal.flux_cals
                + cal.pol_leak_cals
                + cal.pol_angle_cals
            )
        )

    # Remove previous figures
    fnames = glob.glob("{0}/*.png".format(location))
    for fname in fnames:
        os.remove(fname)

    # Generate the plots
    cal.logger.info("Generating plots for manual inspection...")
    plots = []
    for field in fields:

        if fieldtype == "calibrator":
            # Phase vs. Amplitude
            title = "PlotID: {0} Field: {1}".format(len(plots), field)
            plotfile = "{0}/{1}.png".format(location, len(plots))
            casa.plotms(
                vis=cal.vis,
                xaxis="amp",
                yaxis="phase",
                field=field,
                ydatacolumn=datacolumn,
                iteraxis="spw",
                coloraxis="baseline",
                title=title,
                plotfile=plotfile,
                overwrite=True,
                showgui=False,
                exprange="all",
            )
            plots.append(
                {
                    "field": field,
                    "xaxis": "amp",
                    "yaxis": "phase",
                    "avgtime": "",
                    "avgchannel": "",
                }
            )

        # Amplitude vs UV-distance (in wavelength units)
        title = "PlotID: {0} Field: {1}".format(len(plots), field)
        plotfile = "{0}/{1}.png".format(location, len(plots))
        casa.plotms(
            vis=cal.vis,
            xaxis="uvwave",
            yaxis="amp",
            field=field,
            ydatacolumn=datacolumn,
            iteraxis="spw",
            coloraxis="baseline",
            title=title,
            plotfile=plotfile,
            overwrite=True,
            showgui=False,
            exprange="all",
        )
        plots.append(
            {
                "field": field,
                "xaxis": "uvwave",
                "yaxis": "amp",
                "avgtime": "",
                "avgchannel": "",
            }
        )

        # Amplitude vs Time
        title = "PlotID: {0} Field: {1}".format(len(plots), field)
        plotfile = "{0}/{1}.png".format(location, len(plots))
        casa.plotms(
            vis=cal.vis,
            xaxis="time",
            yaxis="amp",
            field=field,
            ydatacolumn=datacolumn,
            iteraxis="spw",
            coloraxis="baseline",
            avgchannel="1e7",
            title=title,
            plotfile=plotfile,
            overwrite=True,
            showgui=False,
            exprange="all",
        )
        plots.append(
            {
                "field": field,
                "xaxis": "time",
                "yaxis": "amp",
                "avgtime": "",
                "avgchannel": "1e7",
            }
        )

        # Amplitude vs Channel
        title = "PlotID: {0} Field: {1}".format(len(plots), field)
        plotfile = "{0}/{1}.png".format(location, len(plots))
        casa.plotms(
            vis=cal.vis,
            xaxis="channel",
            yaxis="amp",
            field=field,
            ydatacolumn=datacolumn,
            iteraxis="spw",
            coloraxis="baseline",
            avgtime="1e7",
            title=title,
            plotfile=plotfile,
            overwrite=True,
            showgui=False,
            exprange="all",
        )
        plots.append(
            {
                "field": field,
                "xaxis": "channel",
                "yaxis": "amp",
                "avgtime": "1e7",
                "avgchannel": "",
            }
        )

        if fieldtype == "calibrator":
            # Phase vs Time
            title = "PlotID: {0} Field: {1}".format(len(plots), field)
            plotfile = "{0}/{1}.png".format(location, len(plots))
            casa.plotms(
                vis=cal.vis,
                xaxis="time",
                yaxis="phase",
                field=field,
                ydatacolumn=datacolumn,
                iteraxis="spw",
                coloraxis="baseline",
                avgchannel="1e7",
                title=title,
                plotfile=plotfile,
                overwrite=True,
                showgui=False,
                exprange="all",
            )
            plots.append(
                {
                    "field": field,
                    "xaxis": "time",
                    "yaxis": "phase",
                    "avgtime": "",
                    "avgchannel": "1e7",
                }
            )

            # Phase vs Channel
            title = "PlotID: {0} Field: {1}".format(len(plots), field)
            plotfile = "{0}/{1}.png".format(location, len(plots))
            casa.plotms(
                vis=cal.vis,
                xaxis="channel",
                yaxis="phase",
                field=field,
                ydatacolumn=datacolumn,
                iteraxis="spw",
                coloraxis="baseline",
                avgtime="1e7",
                title=title,
                plotfile=plotfile,
                overwrite=True,
                showgui=False,
                exprange="all",
            )
            plots.append(
                {
                    "field": field,
                    "xaxis": "channel",
                    "yaxis": "phase",
                    "avgtime": "1e7",
                    "avgchannel": "",
                }
            )
    cal.logger.info("Done.")

    # Generate PDF
    cal.logger.info("Generating PDF...")
    generate_pdf(location)
    cal.logger.info("Done.")

    # Save plot list to a pickle object
    cal.logger.info("Saving plot list to pickle...")
    with open("{0}.pkl".format(location), "w") as fout:
        pickle.dump(plots, fout)
    cal.logger.info("Done.")
